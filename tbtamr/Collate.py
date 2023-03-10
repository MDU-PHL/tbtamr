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

    def _fail_isolate_dict(self,seq_id, step):

        logger.warning(f"Files for the tb-profiler {step} for sample {seq_id} are missing.")
        if step == 'collate':
            logger.critical(f"Will now exit - please check your inputs and try again. You should have files called '{seq_id}/tb-profiler_report.json' at a minimum to proceed with the collation step.")
            raise SystemExit

    def _get_isolate_dict(self, isos):
        # pass
        
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
                    print(isos)
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
        
        Input = namedtuple('Input', 'isolates  exclude_not_reportable')
        to_input = Input(self.isolates, self.exclude_not_reportable) 
        
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
                "Identification (WGS)",
                'Predicted drug resistance',
                "Rifampicin",
                "Rifampicin - Interpretation",
                "Rifampicin - Confidence",
                "Isoniazid",
                "Isoniazid - Interpretation",
                "Isoniazid - Confidence",
                "Pyrazinamide",
                "Pyrazinamide - Interpretation",
                "Pyrazinamide - Confidence",
                "Ethambutol",
                "Ethambutol - Interpretation",
                "Ethambutol - Confidence",
                "Moxifloxacin",
                "Moxifloxacin - Interpretation",
                "Moxifloxacin - Confidence",
                "Amikacin",
                "Amikacin - Interpretation",
                "Amikacin - Confidence",
                "Cycloserine",
                "Cycloserine - Interpretation",
                "Cycloserine - Confidence",
                "Ethionamide",
                "Ethionamide - Interpretation",
                "Ethionamide - Confidence",
                "Kanamycin",
                "Kanamycin - Interpretation",
                "Kanamycin - Confidence",
                "Streptomycin",
                "Streptomycin - Interpretation",
                "Streptomycin - Confidence",
                "Capreomycin",
                "Capreomycin - Interpretation",
                "Capreomycin - Confidence",
                "Para-aminosalicylic acid",
                "Para-aminosalicylic acid - Interpretation",
                "Para-aminosalicylic acid - Confidence",
                "Clofazimine",
                "Clofazimine - Interpretation",
                "Clofazimine - Confidence",
                "Delamanid",
                "Delamanid - Interpretation",
                "Delamanid - Confidence",
                "Bedaquiline",
                "Bedaquiline - Interpretation",
                "Bedaquiline - Confidence",
                "Linezolid",
                "Linezolid - Interpretation",
                "Linezolid - Confidence",
                "Database version"
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

        try:
            with open(f"{path}", 'r') as f:
                results = json.load(f)
                return results
        except Exception as err:
            if for_appending:
                return True
            else:
                logger.critical(f"There seems to have been an error opening {path}. The following error was reported {err}")
    
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

    def _get_interpret(self, drug, mut):

        interp = {'resistance': 'Resistant',
                   'low_level': 'Low-level resistant',
                   'not_reportable': 'Not-reportable',
                   'combo-resistance': 'Resistant only in combination'
                    }
        k = f"{drug}_{mut}"
        
        if 'No mechanism identified' not in mut:
            print(k)
            print(interp[self.db[k]['Confers']])
            return interp[self.db[k]['Confers']]
            
        else:
            return 'Susceptible'
    
    def _get_confidence(self, drug, mut):
        
        k = f"{drug}_{mut}"
        if 'No mechanism identified' not in mut:
            return self.db[k]['Confidence_tbtamr']
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


    def _inference_of_drugs(self, res, drug):
        
        result = []
        inds = [i.strip() for i in res.split(',')]
        for i in inds:
            if i in ["-", "*-","-*"]:
                mut = 'No mechanism identified'
            else:
                mut = self._check_mut(mut = i.strip('*'), drug = drug)
            _d = {  'mutation':mut, 
                    'confidence':self._get_confidence(drug = drug, mut = mut),
                    'interpretation':self._get_interpret(drug = drug, mut = mut)}
            result.append(_d)
        
        return result

    def _infer_drugs(self, tbp_result, seq_id):

        _d = {'Seq_ID':seq_id}
        
        for drug in self.drugs:
            logger.info(f"Checking {drug}")
            _d[self.drugs[drug]] = self._inference_of_drugs(res = tbp_result[seq_id][drug], drug = drug)
        
        return _d
    
    def _infer_dr_profile(self, res):
        
        # print(res)
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

        levs = ['Low-level resistant','Resistant'] if self.exclude_not_reportable else ['Low-level resistant','Resistant','Not-reportable','Resistant only in combination']

        for drug in self.drugs:
            dr = res[self.drugs[drug]]
            
            data = []
            for d in dr:
                
                mt = d['mutation']
                interp = d['interpretation'] if d['interpretation'] in levs else 'Susceptible'
                conf = d['confidence'] if d['interpretation'] in levs else ''
                if mt == 'No mechanism identified':
                    data.append(mt)
                else:
                    data.append(f"{mt} ({interp}|{conf})")
            res[self.drugs[drug]] = ';'.join(data)
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
        for isolate in self.isolates:
            logger.info(f"Working on {isolate}")
            for_collate = self._check_output_file(seq_id=isolate, step = 'collate')
            if for_collate:
                tbp_result = self._open_json(path = self.isolates[isolate]['collate'])
                raw_result = self._open_json(path = self.isolates[isolate]['profile'])
                _dict = self._infer_drugs(tbp_result = tbp_result,seq_id=isolate)
                _dict = self._infer_dr_profile(res = _dict)
                _dict['Species'] = self._species(res = tbp_result, seq_id=isolate)
                _dict['Phylogenetic lineage'] = self._lineage(res = tbp_result, seq_id=isolate)
                _dict['Database version'] = self._db_version(res = raw_result)
                _dict['Median genome coverage'] = self._get_qc_feature(seq_id = isolate, res =tbp_result, val = 'median_coverage')
                _dict['Percentage reads mapped'] = self._get_qc_feature(seq_id = isolate,res = tbp_result, val = 'pct_reads_mapped')
                logger.info(f"Saving results for {isolate}.")
                # saving 'raw' result as json file
                self._save_json(_data = [_dict], _path = f"{isolate}/tbtamr")
                # saving a single sample report - in csv format
                self._single_report(_data = _dict)
                # turn results into a format that is for generic use - larger tabular structure.
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
            mut = ';'.join([i['mutation'] for i in _dict[self.drugs[drug]] if i['interpretation'] in reportable])
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

    def _mdu_infer(self, res):
        
        # print(res)
        

        for drug in self.drugs:
            mut,interpretation,confidence = self._get_infer_conf(_dict = res, drug = drug)
            res[self.drugs[drug]] = mut
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
                # print(_data)
                res.append(_data)
        

        self._save_csv(_data = res, _path = f"{self.sop}_{self.runid}", _type = 'mdu')
        
