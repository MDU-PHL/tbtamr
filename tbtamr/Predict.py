import sys,gzip,pandas,pathlib,json
from CustomLog import logger


import warnings;
warnings.simplefilter(action='ignore', category=FutureWarning)


class PredictAmr(object):

    def __init__(self, 
                 variants,
                 catalog,
                 config,
                 interpretation_rules,
                 classification_rules,
                 seq_id, 
                 vcf,
                 ref,
                 barcode,
                 cascade,
                 call_lineage):
        self.variants = variants
        self.catalog = catalog
        self.config = self.get_config(pth = config)
        self.interpretation_rules = interpretation_rules
        self.seq_id = seq_id
        self.vcf = vcf
        self.ref = ref
        self.barcode = barcode
        self.classification_rules = classification_rules
        self.cascade = cascade
        self.call_lineage = call_lineage
        self.cols = [
            'seq_id',
            'species',
            'main_lineage',
            'predicted drug resistance'
        ] if self.call_lineage else [
            'seq_id',
            'predicted drug resistance'
        ]

    def check_file(self, pth) -> bool:

        if pathlib.Path(pth).exists():
            logger.info(f"{pth} exists.")
        else:
            logger.critical(f"{pth} does not exists or can't be accessed. Please check your inputs and try again.")
            raise SystemExit

        return True
        
    def collect_af(self, variants, var):

        for v in variants:
            if v['variant'] == var:
                return v['af']
        return 0
    
    def collect_resistance_mechs(self,catalog,variants) -> pandas.DataFrame:

        vars = [var['variant'] for var in variants]
        
        mechs = catalog[catalog[self.config['variant_col']].isin(vars)]
        mechs['af'] = mechs[self.config['variant_col']].apply(lambda x:self.collect_af(variants=variants, var = x))
        
        
        return mechs

    def get_highest_conf(self,conf) -> str:

        tp = ""
        if conf != []:
            cfs = sorted([self.config['confidence_levels'][c] for c in conf if c in self.config['confidence_levels']])
            if cfs != []:
                for c in self.config['confidence_levels']:
                    if self.config['confidence_levels'][c] == cfs[0]:
                        tp = self.config['confidence_key'][c]
            else:
                tp = ';'.join(list(set([ self.config['confidence_key'][i] for i in conf])))
        return tp 

    def get_rules_for_dr(self,dr, rules) -> pandas.DataFrame:

        rl_tbl = rules[rules['drug']== dr]
        rl_tbl = rl_tbl.fillna('')
        
        return rl_tbl

    def check_shape(self, rule) -> bool:

        chck = False if rule == "" else True

        return chck

    def extract_mutations(self,dr, result) -> list:
        
        mt = []
        if f"{dr} - mechanisms" in result:
            mt = [m.split()[0] for m in result[f"{dr} - mechanisms"].split(';') if m != 'No reportable mechanims' and self.check_conf_reporting(val = m.split()[-1].strip('()')) ]
            
        return mt

    def update_result(self, result, dr, rule)-> dict:

        if f"{dr} - interpretation" in result:
            new_interp = rule[1]['interpretation']
            result[f"{dr} - interpretation"] = new_interp
            result[f"{dr} - override"] = rule[1]['description'] # add description for tracking purposes.

        return result

    def construct_rule(self, row):
        
        if isinstance(row, tuple):
            d = row[1].to_dict()
        else:
            d = pandas.DataFrame.to_dict(row, orient='records')[0]
        
        rll = []
        for i in range(int(d['number_conditions'])):
            n = i + 1
            values = d[f'values_{n}'].split(';')
            column= d[f'column_{n}']
            # the comparator to check values
            comparator = d[f'comparator_{n}']
            target = d[f'target_{n}']
            # construct the rle
            rle = f"`{column}` {comparator} {values}" 
            
            rll.append(rle)
        
        rule = f"{d['join']}".join(rll) 
        
        return rule

    def apply_rule_override(self, dr, mechs, rules, result) -> dict:
        # get default rule for the drug       
        # logger.info(f"Checking for override rules")

        rl_tbl = self.get_rules_for_dr(dr= dr, rules = rules)
        # print(rl_tbl)
        if not rl_tbl.empty:
            # get all mechs identified in genes associated with dr
            mt = self.extract_mutations(dr = dr, result = result)
            tbl = mechs[mechs[self.config['drug_name_col']].str.lower() == dr]
            
            # check shape
            for row in rl_tbl.iterrows():
                # print(row)
                tbl_to_check = pandas.DataFrame()
                if row[1]['rule_type'] == 'override_simple':
                    # print(f'simple override for {dr}')
                    tbl_to_check = tbl[tbl[self.config['variant_col']].isin(mt)]
                    # if dr == 'rifampicin':
                        # print(tbl_to_check)
                else:
                    tbl_to_check = tbl 
                # print(tbl_to_check)
                rle = self.construct_rule(row = row) 
                # print(rle) 
                # set up date to False by default - leave existing result
                update = False
                # if there is a shape/length criteria to the rule
                if self.check_shape(row[1]['shape']):
                    # logger.info(f"Will check shape of df")
                    shape_rule = f"len(mt) {row[1]['shape']}"
                    if eval(shape_rule): # if the shape/size criteria is met then apply the rest of the rule
                        # print(tbl_to_check)
                        tbl_to_check = tbl_to_check.query(rle) 
                        update = True if not tbl_to_check.empty else False       
                else: # if no size/shape criteria - just apply the rule
                    tbl_to_check = tbl_to_check.query(rle)
                    update = True if not tbl_to_check.empty else False
                # update or not the interpretation.
                result = self.update_result(result = result, dr = dr, rule = row) if update else result 

        return result

    
    def calculate_confidence(self, query_result, m) -> str:

        cf = query_result[query_result[self.config['variant_col']] == m][self.config['confidence_column']].values[0]
        return f"{m} ({self.config['confidence_key'][cf]})"

    def get_resistance_level(self, interp) -> str:

        val = "Susceptible"
        if len(set(interp)) == 1:
            val = interp[0]
        elif len(set(interp)) >1:
            wghts = sorted([self.config['resistance_levels'][i] for i in interp])
            for i in self.config['resistance_levels']:
                if self.config['resistance_levels'][i] == wghts[0]:
                    val = self.config['resistance_levels'][i]
        return val


    def apply_rule_default(self,dr, mechs, rules, result) -> dict: 
        
        # get default rule for the drug
        rl_tbl = self.get_rules_for_dr(dr= dr, rules = rules)
        # get all mechs identified in genes associated with dr
        tbl = mechs[mechs[self.config['drug_name_col']].str.lower() == dr]
        
        muts = []
        conf = []
        interp = []
        bck = []
        # print(dr)
        for row in rl_tbl.iterrows():
        # use rules to identify the appropriate values.
            rle = self.construct_rule(row = row)
        # query the table for the rule
            query_result = tbl.query(rle)
            # extract list of mutations which fit the criteria
            mts = list(query_result[self.config['variant_col']].unique())
            if mts != []:
                mts = [self.calculate_confidence(query_result=query_result, m = m) for m in mts]
            muts.extend(mts)
            conf.extend(list(query_result[self.config['confidence_column']].unique())) 
            # extract interpretation from rules
            if not query_result.empty and row[1]['interpretation'] and row[1]['interpretation'] in self.config['resistance_levels']:
                interp.append(row[1]['interpretation'])
            elif not query_result.empty and row[1]['interpretation'] and not row[1]['interpretation'] in self.config['resistance_levels']:
                # print(row[1]['interpretation'])
                bck.append(row[1]['interpretation'])
        # populate dictionary
        interp = bck if interp == [] else interp # if there are any interpretations that are not in the reported field then they will be sueprseded by the reportable values. If no reportable values are returned then the non reportable may be included - this prevents reporting two or more potentially conflicting results
        hc = self.get_highest_conf(conf = conf)
        ntrp = self.get_resistance_level(interp = interp)
        result[f"{dr} - mechanisms"] = ';'.join(muts) if muts != [] else "No mechanisms identified"
        result[f"{dr} - confidence"] = hc
        result[f"{dr} - interpretation"] = ntrp

        return result
    
    def get_resistance_profile(self, result) -> dict:
        drs = {"first-line":[],"other":[]}
        alldrs = []
        for dt in drs:
            for dr in self.config["drugs_to_report"][dt]:
                if result[f"{dr} - interpretation"] in self.config["resistance_levels"]:
                    drs[dt].append(dr)
                    alldrs.append(dr)
        
        return drs,alldrs
    def get_dlm(self, cond) -> str:
        
        dl = ("","")
        if "&" in cond:
            dl = ("&", " and ")
        elif "|" in cond:
            dl = ("|"," or ")
        return dl
    
    def construct_classification(self, drs, cmprtr, alldrs, dl) -> str:

        drs = drs.split(dl[0]) if dl[0] != "" else [drs]
        req_cond = []
        for rq in drs:
            if rq != "":
                c = f"'{rq}' {cmprtr} {alldrs}"
                req_cond.append(c)
        jn = dl[1]

        return jn.join(req_cond) 

    def get_classification_rule(self, row)-> str:

        lngth = row[1]['shape']
        drg_cls = row[1]['drug_class_condition']
        rq_dl = self.get_dlm(cond = row[1]['required_condition'])
        x_dl = self.get_dlm(cond = row[1]['exlusionary_condition'])
        frst_cond = ""
        if drg_cls != "":
            frst_cond = f"len(drs['{drg_cls}']) {lngth} and"
        rq_rl = self.construct_classification(drs=row[1]['required_condition'], dl = rq_dl, alldrs=f"alldrs", cmprtr=row[1]['comparator'])
        ex_rl = self.construct_classification(drs = row[1]['exlusionary_condition'],cmprtr=row[1]['exclusion_comparator'],alldrs=f"alldrs",dl = x_dl)
        rl = f"{frst_cond} {rq_rl} and ({ex_rl})" if row[1]['exlusionary_condition'] != "" else f"{frst_cond} {rq_rl}"
        
        return rl
        

    def classification(self, rules, result) -> dict:

        drs,alldrs = self.get_resistance_profile(result = result)
        result['predicted drug resistance'] ="No first-line drug resistance"
        for row in rules.iterrows():
            rle = self.get_classification_rule(row = row)
            if eval(rle):
                result['predicted drug resistance'] = row[1]["classification"]
                break

        return result
        
    def compare_mechs_rules(self,interpretation_rules, classification_rules, mechs, result) -> dict:
        
        # print(rules[rules['rule_type'] != 'default'])
        for dr in self.config["drugs_to_infer"]:
            result = self.apply_rule_default(dr = dr.lower(), mechs=mechs, rules=interpretation_rules[interpretation_rules['rule_type'] == 'default'], result = result)
            result = self.apply_rule_override(dr= dr.lower(), mechs=mechs,rules=interpretation_rules[interpretation_rules['rule_type'] != 'default'], result=result)
        result = self.classification(rules = classification_rules, result=result)
        result_df = pandas.DataFrame.from_dict(result, orient= 'index').T
        result_df.to_csv(f"{self.seq_id}/tbtamr_results.csv", index = False)
        return result

    def check_conf_reporting(self, val) -> bool:

        for i in self.config['confidence_key']:
            if self.config['confidence_key'][i] == val and i in self.config['confidence_levels']:
                return True
        return False
    
    def check_for_cascade(self, result, cols) -> bool:

        for c in cols:
            if result[f"{c} - interpretation"] in self.config['resistance_levels']:
                return True
        return False
    
    def generate_drug_cols(self, dr) -> list:

        return [f"{dr} - {i}" for i in ['mechanisms','interpretation','confidence']]

    def cascade_report(self, result, starter_cols):

        # get first level
        
        for dr in self.config['cascade_reporting']['default']:

            starter_cols.extend(self.generate_drug_cols(dr = dr))


        for level in self.config['cascade_reporting']:
            if 'level' in level:
                cols_to_check = self.config['cascade_reporting'][level]["resistance_to_any"]
                if self.check_for_cascade(result=result, cols = cols_to_check):
                    for dr in self.config['cascade_reporting'][level]["report"]:
                        starter_cols.extend(self.generate_drug_cols(dr = dr))
        
        
        return starter_cols
    
    def generate_reporting_df(self, result, output, cols) -> pandas.DataFrame:

        df = pandas.DataFrame.from_dict(result, orient= 'index').T
        df = df[cols]
        rnm = {}
        for c in cols:
            if c == "main_lineage":
                rnm[c] = "Phylogenetic lineage"
            else:
                rnm[c] = c.capitalize()
        df = df.rename(columns=rnm)

        df.to_csv(f"{self.seq_id}/{output}.csv", index = False)        
        
    def make_cascade(self, result) -> bool:

        
        cols = self.cascade_report(result=result, starter_cols = self.cols)
        self.generate_reporting_df(result=result, cols = cols, output="tbtamr_linelist_cascade_report")
    
    
    def make_line_list(self, result, cols) -> bool:

        # wrangle reportable/not reportable
        for dr in self.config['drugs_to_infer']:
            mch = result[f"{dr} - mechanisms"].split(';')
            mchs = []
            # check if mech should be reported - based on conf in string
            for m in mch:
                if m != "No mechanisms identified":
                    c = m.split()[-1].strip('()')
                    if self.check_conf_reporting(val = c):
                        mchs.append(m.split()[0])
                else:
                    mchs.append(m)
            result[f"{dr} - mechanisms"] = ';'.join(mchs)
            # check if conf should be reported
            cf = result[f"{dr} - confidence"].split(';')
            conf = ""
            for c in cf:
                if self.check_conf_reporting(val=c):
                    conf = c
            result[f"{dr} - confidence"] = conf
            # check interpretation
            if result[f"{dr} - interpretation"] not in self.config['resistance_levels']:
                result[f"{dr} - interpretation"] = "Susceptible"
        
        # get cols
        dr2report = self.config['drugs_to_report']
        order = sorted([ o for o in dr2report ])
        
        for o in order:
            c = dr2report[o]
            for dr in c:
                cl = self.generate_drug_cols(dr = dr)
                cols.extend(cl)       
        
        self.generate_reporting_df(result=result, cols = cols, output="tbtamr_linelist_report")
            
        return True

    def get_config(self, pth) -> dict:

        if self.check_file(pth= pth):
            with open(pth, 'r') as j:

                return json.load(j)

    def get_catalog(self, catalog) -> pandas.DataFrame:

        return pandas.read_csv(catalog, dtype= str)

    def get_rules(self, rules) -> pandas.DataFrame:

        rls = pandas.read_csv(rules)
        rls = rls.fillna('')
        return rls
    
    def initiate_results(self) -> dict:

        return {"seq_id": self.seq_id}

    def species(self,lineage) -> str:
        
        species = "Mycobacterium tuberculosis"
        if 'La1' in lineage['main_lineage'] and 'BCG' not in lineage['sub_lineage']:
            species = f"{species} var bovis"
        elif 'La1' in lineage['main_lineage'] and 'BCG' in lineage['sub_lineage']:
            species = f"{species} var bovis BCG"
        elif 'La3' in lineage['main_lineage']:
            species = f"{species} var orygis"
        elif 'La2' in lineage['main_lineage']:
            species = f"{species} var caprae"
        elif lineage['main_lineage'] == '':
            species = "Unable to confirm species"
        
        return species

    def get_main_lineage(self, lins) -> list:

        return list(set([i.split('.')[0] for i in lins]) )

    def wrangle_lineages(self,bca) -> dict:

        lins = sorted([i.id for i in bca])
        main = self.get_main_lineage(lins = lins)
        sl = []
        for m in main:
            sub = ""
            for i in range(1,len(lins)):
                if lins[i-1] in lins[i] and m in lins[i]:
                    sub = lins[i]
            sl.append(sub)

        result =  {
            'main_lineage':';'.join(main), 
            'sub_lineage':';'.join(sl),
        }
        
        result['species'] = self.species(lineage=result)
        return result

    def extract_lineage(self,vcf, brcd, ref, result) -> dict:

        if self.call_lineage:
            from pathogenprofiler import barcode, Vcf
            logger.info("First determining lineage.")
            v = Vcf(vcf)
            muts = v.get_bed_gt(brcd,ref)
            bca = barcode(muts, brcd)
            result = self.wrangle_lineages(bca=bca, result = result)
        
        return result

    def run_prediction(self) -> None:
        
        if self.check_file(pth = self.catalog) and self.check_file(pth= self.classification_rules) and self.check_file(pth= self.interpretation_rules):
            result = self.initiate_results()
            result = self.extract_lineage( vcf = self.vcf, brcd = self.barcode, ref = self.ref, result = result )
            ctlg = self.get_catalog(catalog= self.catalog)
            mechs = self.collect_resistance_mechs(catalog=ctlg, variants=self.variants)
            interpretation_rules = self.get_rules(rules = self.interpretation_rules)
            classification_rules = self.get_rules(rules = self.classification_rules)
            result = self.compare_mechs_rules(interpretation_rules = interpretation_rules, classification_rules=classification_rules,mechs=mechs, result = result)
            self.make_line_list(result = result, cols = self.cols)
            if self.cascade:
                self.make_cascade(result=result)