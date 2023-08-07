#!/usr/bin/env python3
import json

from os import path
from collections import namedtuple
import pandas, pathlib
pandas.options.mode.chained_assignment = None


from tbtamr.CustomLog import logger
from tbtamr.TbTamr import Tbtamr

class Parse(Tbtamr):

    def __init__(self,args):
        super().__init__()
        self.isolates = self._extract_isolates(args.isolates)
        self.exclude_not_reportable = args.exclude_not_reportable
        self.prop_mtb = args.prop_mtb
        self.min_cov = args.min_cov

    def _fail_isolate_dict(self,seq_id, step):

        logger.warning(f"Files for the tb-profiler {step} for sample {seq_id} are missing.")
        if step == 'collate':
            logger.critical(f"Will now exit - please check your inputs and try again. You should have files called '{seq_id}/tb-profiler_report.json' at a minimum to proceed with the collation step.")
            raise SystemExit

    def _get_isolate_dict(self, isos):
        
        isolates = {}
        for i in isos:
            logger.info(f"Checking if files for {i} are present.")
            isolates[i] = {}
            pr = self._check_output_file(seq_id=i, step='profile')
            isolates[i]['profile'] = pr if pr else self._fail_isolate_dict(seq_id = i,step = 'profile')
            cr = self._check_output_file(seq_id=i, step='collate')
            isolates[i]['collate'] = cr if cr else self._fail_isolate_dict(seq_id = i,step = 'collate')
        
        return isolates

    def _extract_isolates(self, _input):

        if self._file_present(name = _input) and pathlib.Path(_input).is_dir():
            isolates = self._get_isolate_dict(isos = [_input])
            return isolates
        
        elif self._file_present(name = _input) and not pathlib.Path(_input).is_dir():
            try:
                with open(_input, 'r') as f:
                    isos = f.read().strip().split('\n')
                    if isinstance(isos, list) and len(isos)>1:
                        isolates = self._get_isolate_dict(isos = isos)
                        return isolates
                    else:
                        logger.critical(f"It seems that your input file is not configured properly. Isolates should be listed with a new isolate on each line. Please try again.")
                        raise SystemExit
            except Exception as err:
                logger.critical(f"Was unable to open {_input} and extract isolates. The following error was reported {err}")
                raise SystemExit
    
    def extract_inputs(self):
        
        Input = namedtuple('Input', ['isolates'  ,'exclude_not_reportable', 'prop_mtb','min_cov'])
        to_input = Input(self.isolates, self.exclude_not_reportable, self.prop_mtb, self.min_cov) 
        
        return to_input

class Inferrence(Tbtamr):

    """
    a class for collation of tbprofiler results.
    """
    
    def __init__(self,args):
        super().__init__()
        self.isolates = args.isolates
        self.db_path = f"{pathlib.Path(__file__).parent.parent /'tbtamr' /'db' / 'tbtamr_db_latest.json'}"
        self.db = self._get_db(path = self.db_path)
        self.drugs = self._get_drugs()
        self.low_level = self._get_low_level()
        self.not_reportable = self._get_not_reportable()
        self.exclude_not_reportable = args.exclude_not_reportable
        self.prop_mtb = args.prop_mtb
        self.min_cov = args.min_cov
        
    def _get_db(self,path):

        try:
            with open(path, 'r') as j:
                return json.load(j)
        except Exception as err:
            logger.critical(f"There seems to have been an error opening {path}. The following error was reported {err}")


    def _get_low_level(self):
        
        low_level = []
        for record in self.db:
            if self.db[record]['Confers'] == 'low_level':
                low_level.append(record)
        
        return low_level

    def _get_not_reportable(self):

        not_reportable = []
        for record in self.db:
            if self.db[record]['Confers'] == 'not_reportable':
                not_reportable.append(record)
        
        return not_reportable

    def _get_header(self,_type):

        if _type == 'mdu':
            header =  [
                "Seq_ID",
                "Species",
                "Phylogenetic lineage",
                'Predicted drug resistance',
                "Rifampicin - ResMech",
                "Rifampicin - Interpretation",
                "Rifampicin - Confidence",
                "Isoniazid - ResMech",
                "Isoniazid - Interpretation",
                "Isoniazid - Confidence",
                "Pyrazinamide - ResMech",
                "Pyrazinamide - Interpretation",
                "Pyrazinamide - Confidence",
                "Ethambutol - ResMech",
                "Ethambutol - Interpretation",
                "Ethambutol - Confidence",
                "Moxifloxacin - ResMech",
                "Moxifloxacin - Interpretation",
                "Moxifloxacin - Confidence",
                "Amikacin - ResMech",
                "Amikacin - Interpretation",
                "Amikacin - Confidence",
                "Cycloserine - ResMech",
                "Cycloserine - Interpretation",
                "Cycloserine - Confidence",
                "Ethionamide - ResMech",
                "Ethionamide - Interpretation",
                "Ethionamide - Confidence",
                "Para-aminosalicylic acid - ResMech",
                "Para-aminosalicylic acid - Interpretation",
                "Para-aminosalicylic acid - Confidence",
                "Kanamycin - ResMech",
                "Kanamycin - Interpretation",
                "Kanamycin - Confidence",
                "Streptomycin - ResMech",
                "Streptomycin - Interpretation",
                "Streptomycin - Confidence",
                "Capreomycin - ResMech",
                "Capreomycin - Interpretation",
                "Capreomycin - Confidence",
                "Clofazimine - ResMech",
                "Clofazimine - Interpretation",
                "Clofazimine - Confidence",
                "Delamanid - ResMech",
                "Delamanid - Interpretation",
                "Delamanid - Confidence",
                "Bedaquiline - ResMech",
                "Bedaquiline - Interpretation",
                "Bedaquiline - Confidence",
                "Linezolid - ResMech",
                "Linezolid - Interpretation",
                "Linezolid - Confidence",
                "Quality",
                "Database version",
            ]
        else:
            header =  [
                "Seq_ID",
                "Species",
                "Phylogenetic lineage",
                'Predicted drug resistance',
                "Rifampicin",
                "Isoniazid",
                "Pyrazinamide",
                "Ethambutol",
                "Moxifloxacin",
                "Amikacin",
                "Cycloserine",
                "Ethionamide",
                "Kanamycin",
                "Streptomycin",
                "Capreomycin",
                "Para-aminosalicylic acid",
                "Clofazimine",
                "Delamanid",
                "Bedaquiline",
                "Linezolid",
                "Median genome coverage",
                "Percentage reads mapped",
                "Quality",
                "Database version"
            ]

        return header


    def _get_drugs(self):

        drugs = {
            "rifampicin":"Rifampicin",
            "isoniazid":"Isoniazid",
            "pyrazinamide":"Pyrazinamide",
            "ethambutol":"Ethambutol",
            "moxifloxacin":"Moxifloxacin",
            "amikacin":"Amikacin",
            "cycloserine":"Cycloserine",
            "ethionamide":"Ethionamide",
            "para-aminosalicylic_acid":"Para-aminosalicylic acid",
            "clofazimine":"Clofazimine",
            "delamanid":"Delamanid",
            "bedaquiline":"Bedaquiline",
            "linezolid":"Linezolid",
            "kanamycin":"Kanamycin",
            "streptomycin":"Streptomycin",
            "capreomycin":"Capreomycin"
        }

        return drugs
   
    
    def _open_json(self, path, for_appending = False):
        # print(path)
        try:
            with open(f"{path}", 'r') as f:
                results = json.load(f)
                return results
        except Exception as err:
            if for_appending:
                return True
            else:
                logger.critical(f"There seems to have been an error opening {path}. The following error was reported {err}")
                raise SystemExit
            
    def _save_json(self, _data, _path):

        
        logger.info(f"Checking {_path}.json already exists.")
        _existing_data = self._open_json(path = f"{_path}.json",for_appending=True)
        if isinstance(_existing_data, dict):
            _to_save = _existing_data.update(_data)
        else:
            _to_save = _data
            
        logger.info(f"Saving file: {_path}.json")
        with open(f"{_path}.json", 'w') as f:
            json.dump(_to_save,f, indent = 4)

    def _save_csv(self, _data, _path, _type):
        
        logger.info(f"Saving file: {_path}.csv")

        header = self._get_header(_type = _type)

        if isinstance(_data, dict):       
            df = pandas.DataFrame(_data, index = [0])
        else:
            df = pandas.DataFrame(_data)
        
        df = df[header]
        
        df.to_csv(f"{_path}.csv", index = False)

    def _save_data(self, _data, prefix):

        self._save_csv(_data,prefix, 'general')
        self._save_json(_data,prefix)

    def _check_quality(self, cov, perc):
        
        if cov >= self.min_cov and perc >= self.prop_mtb:

            return 'Pass QC'
        elif perc < self.prop_mtb:
            return 'Failed: < 80 % _M. tuberculosis_ reads in sample'
        elif  cov < self.min_cov:
            return 'Failed: < 40x aligned coverage to reference genome'

    def _get_interpret(self, drug, mut, qual):

        if 'Fail' not in qual:
            
            interp = {'resistance': 'Resistant',
                    'low_level': 'Low-level resistant',
                    'not_reportable': 'Not-reportable',
                    'combo-resistance': 'Resistant only in combination'
                        }
            k = f"{drug}_{mut}"
            
            if 'No mechanism identified' not in mut:
                return interp[self.db[k]['Confers']]
                
            else:
                return 'Susceptible'
        else:
            return ''
    
    def _get_confidence(self, drug, mut, qual):
        
        if not 'Fail' in qual:
            logger.info(f"Sequence passed QC.")
            k = f"{drug}_{mut}"
            if 'No mechanism identified' not in mut:
                return self.db[k]['Confidence_tbtamr']
            else:
                return ''
        else:
            return ''

    def _check_mut(self, drug, mut):

        if f"{drug}_{mut}" in self.db:
            return mut
        else:
            gene = mut.split('_')[0]
            if gene in ['ethA','katG','thyA','pncA','gid'] and 'del' in mut:
                return f"{gene}_large_deletion"
            else:
                logger.critical(f"There is a problem with the mutation reported {mut} - it is not in the database of expected mutations. This will not be reported. Please leave an issue on github")
                return "No mechanism identified*"


    def _inference_of_drugs(self, res, drug, qual):
        
        result = []
        inds = [i.strip() for i in res.split(',')]
        for i in inds:
            if i in ["-", "*-","-*"]:
                mut = 'No mechanism identified'
            else:
                mut = self._check_mut(mut = i.strip('*'), drug = drug)
            _d = {  'mutation':mut if 'Fail' not in qual else qual, 
                    'confidence':self._get_confidence(drug = drug, mut = mut, qual = qual),
                    'interpretation':self._get_interpret(drug = drug, mut = mut, qual = qual)}
            # print(_d)
            result.append(_d)
        
        return result

    def _infer_drugs(self, tbp_result, seq_id, qual):

        _d = {'Seq_ID':seq_id}
        
        for drug in self.drugs:
            logger.info(f"Checking {drug}")
            _d[self.drugs[drug]] = self._inference_of_drugs(res = tbp_result[seq_id][drug], drug = drug, qual = qual)
        
        return _d
    
    def _infer_dr_profile(self, res, qual):
        
        if 'Fail' not in qual:
            fline_b = ['pyrazinamide','ethambutol']
            flq = ['ofloxacin','moxifloxacin','levofloxacin']
            other_groupA = ["bedaquiline","linezolid","delamanid"]
            score = 0
            other_score = False
            flq_score = False
            inh = 3
            rif = 8
            fline_score = 1
            rf = ''
            
            for drug in self.drugs:
                # print(drug)
                # print(res[self.drugs[drug]])
                if 'No mechanism identified' not in res[self.drugs[drug]][0]['mutation']:
                    # get the interpretations from the results
                    interp = [i['interpretation'] for i in res[self.drugs[drug]]]
                    if drug in flq and any(item in interp for item in ['Low-level resistant','Resistant']):
                        flq_score = True
                    if drug in other_groupA and any(item in interp for item in ['Low-level resistant','Resistant']):
                        other_score = True
                    if drug == 'rifampicin' and any(item in interp for item in ['Low-level resistant','Resistant']):
                        rf = ' (RR-TB)'
                        score = score + rif
                    elif drug == 'isoniazid' and any(item in interp for item in ['Low-level resistant','Resistant']):
                        score = score + inh
                    elif drug in fline_b and any(item in interp for item in ['Low-level resistant','Resistant']):
                        score = score + fline_score
            resistance = 'No first-line drug resistance predicted'
            if score >=8 and flq_score == True and other_score == False:
                resistance = 'Pre-Extensive/Extensive drug-resistance predicted'
            elif score >=8 and flq_score == True and other_score == True:
                resistance = 'Pre-Extensive/Extensive drug-resistance predicted'
            elif score in [1,3,8]: # one first line drug
                resistance = f'Mono-resistance predicted{rf}'
            elif score in [2,4,5,9,10]: # more than one first line drug where INH OR RIF can be present
                resistance = f'Poly-resistance predicted{rf}'
            elif score >=11 and flq_score == False:
                resistance = 'Multi-drug resistance predicted (MDR-TB)'
            res['Predicted drug resistance'] = f"{resistance}"
        else:
            res['Predicted drug resistance'] = qual

        return res

    def _species(self,res, seq_id):
        
        species = "Mycobacterium tuberculosis"
        if 'La1' in res[seq_id]['main_lin'] and 'BCG' not in res[seq_id]['sublin']:
            species = f"{species} var bovis"
        elif 'La1' in res[seq_id]['main_lin'] and 'BCG' in res[seq_id]['sublin']:
            species = f"{species} var bovis BCG"
        elif 'La3' in res[seq_id]['main_lin']:
            species = f"{species} var orygis"
        elif 'La2' in res[seq_id]['main_lin']:
            species = f"{species} var caprae"
        elif res[seq_id]['main_lin'] == '':
            species = "Not likely M. tuberculosis"
        
        return species

    def _lineage(self,res, seq_id):
        
        if res[seq_id]['main_lin'] == '':
            return "Not typable"
        else:
            return res[seq_id]['main_lin'].replace('lineage', 'Lineage ')

    def _db_version(self, res):
        try:
            return f"{res['db_version']['name']}_{res['db_version']['commit']}"
        except:
            return f"No database version available"
    
    def _get_qc_feature(self,seq_id, res, val):

        return res[seq_id][val]

    def _wrangle_json(self,res):
        # print(res)
        levs = ['Low-level resistant','Resistant'] if self.exclude_not_reportable else ['Low-level resistant','Resistant','Not-reportable','Resistant only in combination']
        # print(res)
        for drug in self.drugs:
            dr = res[self.drugs[drug]]
            
            data = set()
            for d in dr:
                
                mt = d['mutation']
                interp = d['interpretation'] if d['interpretation'] in levs or d['interpretation'] == ''  else 'Susceptible'
                conf = d['confidence'] if d['interpretation'] in levs else ''
                if mt == 'No mechanism identified' or 'Fail' in mt:
                    data.add(mt)
                else:
                    data.add(f"{mt} ({interp}|{conf})")
            res[self.drugs[drug]] = ';'.join(list(data))
        
        logger.info(f"Saving a wide format csv for {res['Seq_ID']}")
        self._save_csv(_data = res, _path = f"{res['Seq_ID']}/tbtamr", _type = 'general')
        
        return res


    def _single_report(self, _data):

        lines = [
            f"Seq_ID:,{_data['Seq_ID']}",
            f"Species:,{_data['Species']}",
            f"Phylogenetic lineage:,{_data['Phylogenetic lineage']}",
            f"Seq_ID:,{_data['Seq_ID']}",
            f"Predicted drug resistance:,{_data['Predicted drug resistance']}",
            f"",
            f"Drug,Mutation,Interpretation,Confidence",
        ]

        for drug in self.drugs:

            dr = _data[self.drugs[drug]]
            for d in dr:
                lines.append(f"{drug.capitalize()},{d['mutation']},{d['interpretation']},{d['confidence']}")
        
        lines.append("")
        lines.append(f"Database version:,{_data['Database version']}")

        logger.info(f"Saving report for {_data['Seq_ID']}")
        pathlib.Path(f"{_data['Seq_ID']}/tbtamr_report.csv").write_text('\n'.join(lines))

    def infer(self):
        
        logger.info(f"Now inferring resistance profiles")
        results = []
        # print(self.isolates)
        for isolate in self.isolates:
            logger.info(f"Working on {isolate}")
            for_collate = self._check_output_file(seq_id=isolate, step = 'collate')
            if for_collate:
                # _dict = {}
                tbp_result = self._open_json(path = self.isolates[isolate]['collate'])
                raw_result = self._open_json(path = self.isolates[isolate]['profile'])
                med_cov = self._get_qc_feature(seq_id = isolate, res =tbp_result, val = 'median_coverage')
                perc_mtb = self._get_qc_feature(seq_id = isolate,res = tbp_result, val = 'pct_reads_mapped')
                qual = self._check_quality(cov = med_cov,perc = perc_mtb)
                
                # print(_dict)
                _dict = self._infer_drugs(tbp_result = tbp_result,seq_id=isolate, qual = qual)
                _dict = self._infer_dr_profile(res = _dict,qual = qual)
                _dict['Species'] = self._species(res = tbp_result, seq_id=isolate) if qual == 'Pass QC' else qual
                _dict['Phylogenetic lineage'] = self._lineage(res = tbp_result, seq_id=isolate) if qual == 'Pass QC' else qual
                _dict['Database version'] = self._db_version(res = raw_result)
                _dict['Quality'] = qual
                _dict['Median genome coverage'] = med_cov
                _dict['Percentage reads mapped'] = perc_mtb
                logger.info(f"Saving results for {isolate}.")
                # saving 'raw' result as json file
                self._save_json(_data = [_dict], _path = f"{isolate}/tbtamr")
                # saving a single sample report - in csv format
                self._single_report(_data = _dict)
                # turn results into a format that is for generic use - larger tabular structure.
                # print(_dict)
                results.append(self._wrangle_json(res = _dict))
            
        logger.info(f"Saving collated data.")
        self._save_csv(_data = results, _path = "tbtamr", _type = 'general')



class Mdu(Inferrence):

    """
    a class for collation of tbprofiler results.
    """
    
    def __init__(self,args):
        # super().__init__()
        self.jsons = args.json
        self.sop = args.output_name
        self.runid = args.runid
        # self.db_path = f"{pathlib.Path(__file__).parent.parent /'tbtamr' /'db' / 'tbtamr_db_latest.json'}"
        # self.db = self._get_db(path = self.db_path)
        self.drugs = self._get_drugs()
        # self.low_level = self._get_low_level()
        # self.not_reportable = self._get_not_reportable()
    
    def _get_infer_conf(self, _dict, drug):

        reportable = ['Low-level resistant','Resistant']
        if _dict['Species'] == 'Mycobacterium tuberculosis':
            interp = [i['interpretation'] for i in _dict[self.drugs[drug]] if i['interpretation'] in reportable]
            conf = [i['confidence'] for i in _dict[self.drugs[drug]]]
            mut = ';'.join([i['mutation'].replace('_',' ') for i in _dict[self.drugs[drug]] if i['interpretation'] in reportable])
            if interp != []:
                interpretation = sorted(interp)[-1]
                confidence = sorted(conf)[0]
            else:
                interpretation = 'Susceptible'
                confidence = ''
                mut = "No mechanism identified"
        else:
            mut = interpretation = 'Not reportable'
            confidence = ''
        
        return mut,interpretation,confidence
    
    def _get_summary_res(self, res):

        if not res['Species'] in ["Mycobacterium tuberculosis", "Failed: < 40x aligned coverage to reference genome", "Failed: < 80 % _M. tuberculosis_ reads in sample"] :
            return "Not reportable"
        else:
            return res["Predicted drug resistance"]

    def _mdu_infer(self, res):
        
        # print(res)
        

        for drug in self.drugs:
            mut,interpretation,confidence = self._get_infer_conf(_dict = res, drug = drug)
            res[f"{self.drugs[drug]} - ResMech"] = mut
            res[f"{self.drugs[drug]} - Interpretation"] = interpretation
            res[f"{self.drugs[drug]} - Confidence"] = confidence
        
        return res
            

            

    def mduify(self):
        headers = self._get_header(_type = 'mdu')
        res = []
        for j in self.jsons:

            if pathlib.Path(j).exists():
                js = self._open_json(j)
                _data = self._mdu_infer(res = js[0])
                _data["Predicted drug resistance"] = self._get_summary_res(res = _data)
                # print(_data)
                res.append(_data)
        

        self._save_csv(_data = res, _path = f"{self.sop}_{self.runid}", _type = 'mdu')
        
