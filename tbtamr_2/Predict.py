import sys,gzip,pandas,pathlib,json
from CustomLog import logger
from pathogenprofiler import barcode, Vcf

import warnings;
warnings.simplefilter(action='ignore', category=FutureWarning)


class PredictAmr(object):

    def __init__(self, 
                 variants,
                 catalog,
                 config,
                 rules,
                 seq_id, 
                 vcf,
                 ref,
                 barcode):
        self.variants = variants
        self.catalog = catalog
        self.config = self.get_config(pth = config)
        self.rules = rules
        self.seq_id = seq_id
        self.vcf = vcf
        self.ref = ref
        self.barcode = barcode

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
            cfs = sorted([self.config['confidence_levels'][c] for c in conf])
            
            for c in self.config['confidence_levels']:
                if self.config['confidence_levels'][c] == cfs[0]:
                    tp = self.config['confidence_key'][c]

        return tp 

    def get_rules_for_dr(self,dr, rules, rule_type) -> pandas.DataFrame:
        rl_tbl = rules[rules['drug']== dr]
        
        rl_tbl = rl_tbl.fillna('')
        # check that table is not empty
        if rl_tbl.shape[0] != 1 and rule_type == 'default':
            logger.critical(f"You must have one default rule. Please check your inputs and try again.")
            raise SystemExit
        return rl_tbl

    def check_shape(self, rule) -> bool:

        chck = False if rule == "" else True

        return chck

    def extract_mutations(self,dr, result) -> list:
        mt = []
        if f"{dr} - mechanisms" in result:
            mt = [m for m in result[f"{dr} - mechanisms"].split(';') if m != 'No reportable mechanims' ]
            
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
        # get all mechs identified in genes associated with dr
        rl_tbl = self.get_rules_for_dr(dr= dr, rules = rules, rule_type= 'overrride')
        if not rl_tbl.empty:
            # get all mechs identified in genes associated with dr
            mt = self.extract_mutations(dr = dr, result = result)
            tbl = mechs[mechs[self.config['drug_name_col']].str.lower() == dr]
            # check shape
            for row in rl_tbl.iterrows():
                if row[1]['rule_type'] == 'override_simple':
                    tbl_to_check = tbl[tbl[self.config['variant_col']].isin(mt)]
                else:
                    tbl_to_check = tbl
                rle = self.construct_rule(row = row)
                # set up date to False by default - leave existing result
                update = False
                # if there is a shape/length criteria to the rule
                if self.check_shape(row[1]['shape']):
                    # logger.info(f"Will check shape of df")
                    shape_rule = f"len(mt) {row[1]['shape']}"
                    if eval(shape_rule): # if the shape/size criteria is met then apply the rest of the rule
                        tbl_to_check = tbl_to_check.query(rle) 
                        update = True if not tbl_to_check.empty else False
                        
                else: # if no size/shape criteria - just apply the rule
                    tbl_to_check = tbl_to_check.query(rle)
                    update = True if not tbl_to_check.empty else False
                # update or not the interpretation.
                result = self.update_result(result = result, dr = dr, rule = row) if update else result 

        return result

    

    def apply_rule_default(self,dr, mechs, rules, result) -> dict: 
        
        # get default rule for the drug
        rl_tbl = self.get_rules_for_dr(dr= dr, rules = rules, rule_type= 'default')
        # get all mechs identified in genes associated with dr
        tbl = mechs[mechs[self.config['drug_name_col']].str.lower() == dr]
        # use rules to identify the appropriate values.
        rle = self.construct_rule(row = rl_tbl)
        # query the table for the rule
        tbl = tbl.query(rle)
        # extract list of mutations which fit the criteria
        muts = ';'.join(list(tbl[self.config['variant_col']].unique())) if not tbl.empty else 'No reportable mechanims'
        tbl['mutaf'] = tbl[[self.config['variant_col'], 'af']].apply(lambda x : f"{x[0]} ({x[1]})", axis = 1)
        mutaf = ';'.join(list(tbl['mutaf'].unique())) if not tbl.empty else ''
        # calculate highest confidence
        conf = list(tbl[self.config['confidence_column']].unique()) if not tbl.empty else []
        hc = self.get_highest_conf(conf = conf)
        # extract interpretation from rules
        interp = rl_tbl['interpretation'].values[0] if not tbl.empty else 'Susceptible'
        # populate dictionary
        result[f"{dr} - mechanisms"] = muts 
        result[f"{dr} - mechanisms (freq)"] = mutaf 
        result[f"{dr} - confidence"] = hc 
        result[f"{dr} - interpretation"] = interp

        return result

    def compare_mechs_rules(self,rules, mechs, result) -> True:
        
        # print(rules[rules['rule_type'] != 'default'])
        for dr in self.config["drugs_to_infer"]:
            result = self.apply_rule_default(dr = dr.lower(), mechs=mechs, rules=rules[rules['rule_type'] == 'default'], result = result)
            result = self.apply_rule_override(dr= dr.lower(), mechs=mechs,rules=rules[rules['rule_type'] != 'default'], result=result)
        result_df = pandas.DataFrame.from_dict(result, orient= 'index').T
        result_df.to_csv(f"{self.seq_id}/tbtamr_results.csv", index = False)

    def get_config(self, pth):

        if self.check_file(pth= pth):
            with open(pth, 'r') as j:

                return json.load(j)

    def get_catalog(self, catalog) -> pandas.DataFrame:

        return pandas.read_csv(catalog, dtype= str)

    def get_rules(self, rules) -> pandas.DataFrame:

        return pandas.read_csv(rules)
    
    def species(self,lineage):
        
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
            'seq_id': self.seq_id,
            'main_lineage':';'.join(main), 
            'sub_lineage':';'.join(sl),
        }
        
        result['species'] = self.species(lineage=result)
        return result

    def extract_lineage(self,vcf, brcd, ref) -> dict:

        logger.info("First determining lineage.")
        v = Vcf(vcf)
        muts = v.get_bed_gt(brcd,ref)
        bca = barcode(muts, brcd)
        lineage = self.wrangle_lineages(bca)
        
        return lineage

    def run_prediction(self) -> True:
        
        if self.check_file(pth = self.catalog) and self.check_file(pth= self.rules):
            lineage = self.extract_lineage( vcf = self.vcf, brcd = self.barcode, ref = self.ref )
            ctlg = self.get_catalog(catalog= self.catalog)
            mechs = self.collect_resistance_mechs(catalog=ctlg, variants=self.variants)
            rules = self.get_rules(rules = self.rules)

            if self.compare_mechs_rules(rules = rules, mechs=mechs, result = lineage):
                logger.info(f"Prediction for {self.seq_id} has been completed and saved as {self.seq_id}/tbtamr_result_{self.seq_id}.csv")