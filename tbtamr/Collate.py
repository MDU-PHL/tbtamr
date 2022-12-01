#!/usr/bin/env python3
import json
from os import path
import pandas,pathlib
pandas.options.mode.chained_assignment = None

# from pandas.core.algorithms import isin
from tbtamr.CustomLog import logger
from tbtamr.TbTamr import Tbtamr

class Parse(Tbtamr):

    def __init__(self,args):
        super().__init__()
        self.isolates = self._extract_isolates(args.isolates)
    
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

        return self.isolates

class Inferrence(Tbtamr):

    """
    a class for collation of tbprofiler results.
    """
    
    def __init__(self,args):
        super().__init__()
        self.isolates = args
        self.drugs = self._get_drugs()

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
   
    def _get_funcs(self):
        funcs= {
                'rifampicin':self._rifampicin,
                'isoniazid': self._isoniazid,
                'pyrazinamide': self._pyrazinamide,
                'ethambutol':self._ethambutol,
                'moxifloxacin':self._other_drugs,
                'ofloxacin':self._other_drugs,
                'amikacin':self._other_drugs,
                'capreomycin':self._other_drugs,
                'kanamycin':self._other_drugs,
                'ethionamide':self._other_drugs,
                'streptomycin':self._other_drugs,
                'cycloserine':self._other_drugs,
                'para-aminosalicylic_acid': self._other_drugs,
                'clofazimine': self._other_drugs,
                'delamanid':self._other_drugs,
                'bedaquiline':self._other_drugs,
                'linezolid':self._other_drugs,
                'kanamycin':self._other_drugs,
                'streptomycin':self._other_drugs,
                'capreomycin':self._other_drugs
                }

        return funcs

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
            json.dump(_to_save,f)

    def _save_csv(self, _data, _path):
        
        logger.info(f"Saving file: {_path}.csv")

        header = [
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

        
        df = pandas.DataFrame(_data)
        df = df[header]
        
        df.to_csv(f"{_path}.csv", index = False)

    def _save_data(self, _data, prefix):

        self._save_csv(_data,prefix)
        self._save_json(_data,prefix)


    def _rifampicin(self, res):
        rif_to_remove = [
            'rpoB_p.Ala286Val',
            'rpoB_p.Arg552Cys',
            'rpoB_p.Glu423Ala',
            'rpoB_p.Glu460Gly',
            'rpoB_p.Glu592Asp',
            'rpoB_p.Glu761Asp',
            'rpoB_p.Gly981Asp',
            'rpoB_p.Ile480Thr',
            'rpoB_p.Ile480Val',
            'rpoB_p.Leu490Val',
            'rpoB_p.Phe424Leu',
            'rpoB_p.Phe424Val',
            'rpoB_p.Pro454Arg',
            'rpoB_p.Pro454His',
            'rpoB_p.Pro454Leu',
            'rpoB_p.Pro454Ser',
            'rpoB_p.Pro483Leu',
            'rpoB_p.Ser493Leu',
            'rpoB_p.Thr400Ala',
            'rpoB_p.Gly332Arg'
            ]
        rif_borderline = ["rpoB_p.Leu430Pro","rpoB_p.His445Asn","rpoB_p.Asp435Tyr",
                      "rpoB_p.Leu452Pro","p.Ile491Phe"]
#     print(x)
        if res == '-':
            return 'No mechanism identified'
        elif res in ["*-","-*"]:
            return 'No mechanism identified^'
        else:
            inds = [i.strip() for i in res.split(',')]
#         print(inds)
        if len(inds) == 1 and inds[0] in rif_to_remove:
            return 'No mechanism identified'
        elif len(inds) == 1 and inds[0] in rif_borderline:
            return f'{inds[0]}*'
        else:
            return f";".join(inds)

    def _isoniazid(self,res):
        if res == '-':
            return 'No mechanism identified'
        elif res in ["*-","-*"]:
            return 'No mechanism identified^'
        else:
            inds = [i.strip() for i in res.split(',')]
        if len(inds) == 1 and 'ahpC' in inds[0]:
            return 'No mechanism identified'
        elif len(inds) == 1 and 'fabG' in inds[0]:
            return f'{inds[0]}*' if inds[0] == 'fabG1_c.-15C>T' else 'No mechanism identified'
        else:
            return f";".join(inds)

    def _emb_comp(self, res):
    
        gln497_hc = ['embB_p.Gln497Arg','embB_p.Gln497Lys']
        leu370_us = ['embB_p.Leu370Arg']
        embA = ['embA_c.-16C>T']
        to_remove = []
        for r in res:
            if 'embB_p.Gln497' in r and r not in gln497_hc:
                to_remove.append(r)
            elif r in leu370_us:
                to_remove.append(r)
            elif r in embA:
                to_remove.append(r)
        if len(set(res).difference(to_remove))== 0:
            return 'No mechanism identified'
        else:
            return ';'.join(list(set(res).difference(to_remove)))

    def _ethambutol(self,res):
             
        if res == '-':
            return 'No mechanism identified'
        elif res in ["*-","-*"]:
            return 'No mechanism identified^'
        else:
            inds = [i.strip() for i in res.split(',')]
            return self._emb_comp(res = inds)
    
    def _remove_who(self,res):

        not_reportable = ['PPE35','clpC1','Rv3236c', 'Rv1258c','panD']
        to_remove = []
        for i in res:
            for r in not_reportable:
                if r in i:
                    to_remove.append(i)
        
        if len(set(res).difference(to_remove))== 0:
            return 'No mechanism identified*'
        else:    
            return f";".join(list(set(res).difference(to_remove)))

    def _pyrazinamide(self,res):
        
        if res == '-':
            return 'No mechanism identified'
        elif res in ["*-","-*"]:
            return 'No mechanism identified^'
        inds = [i.strip() for i in res.split(',')]
        
        return self._remove_who(inds)

    def _other_drugs(self, res):
        if res == '-':
            return 'No mechanism identified'
        elif res in ["*-","-*"]:
            return 'No mechanism identified^'
        return f";".join(res.split(','))
    
    def _infer_drugs(self, tbp_result, seq_id):

        _d = {'Seq_ID':seq_id}
        funcs = self._get_funcs()
        for drug in self.drugs:
            logger.info(f"Checking {drug}")
            _d[self.drugs[drug]] = funcs[drug](tbp_result[seq_id][drug])
        
        return _d
    
    def _infer_dr_profile(self, res):
        
        # fline_a = ['rifampicin','isoniazid']
        fline_b = ['pyrazinamide','ethambutol']
        flq = ['ofloxacin','moxifloxacin','levofloxacin']
        other_groupA = ["bedaquiline",	"linezolid","delamanid"]
        score = 0
        other_score = False
        flq_score = False
        inh = 3
        rif = 8
        fline_score = 1
        rf = ''
        
        suff = ""
        for drug in self.drugs:
            if 'No mechanism identified' not in res[self.drugs[drug]]:
                if drug in flq:
                    flq_score = True
                if drug in other_groupA:
                    other_score = True
                if drug == 'rifampicin':
                    rf = ' (RR-TB)'
                    score = score + rif
                elif drug == 'isoniazid':
                    score = score + inh
                elif drug in fline_b:
                    score = score + fline_score
                
            if "*" in res[self.drugs[drug]] : #if there are any intermediate resistance mechanims identified add a *
                suff = f"{suff}*"
            elif "^" in res[self.drugs[drug]]: # if there are any missing
                suff = f"{suff}^"
        print(f"this sample has a score of {score}, flq is {flq_score} and other is {other_score}") 
        resistance = 'No drug resistance predicted'
        if score >=8 and flq_score == True and other_score == False:
            resistance = 'Pre-Extensive/Extensive drug-resistance predicted'
        elif score >=8 and flq_score == True and other_score == True:
            resistance = 'Pre-Extensive/Extensive drug-resistance predicted'
        elif score in [1,3,8]: # one first line drug
            resistance = f'Mono-resistance predicted{rf}'
        elif score in [2,4,5,9,10]: # more than one first line drug where INH OR RIF can be present
            resistance = f'Poly-resistance predicted{rf}'
        elif score >=11 and flq_score == False and other_score == False:
            resistance = 'Multi-drug resistance predicted (MDR-TB)'
        
                       

        res['Predicted drug resistance'] = f"{resistance}{suff}"
        
        return res

    def _species(self,res, seq_id):
        # pass
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


    def infer(self):
        
        logger.info(f"Now inferring resistance profiles")
        results = []
        for isolate in self.isolates:
            logger.info(f"Working on {isolate}")
            for_collate = self._check_output_file(seq_id=isolate, step = 'collate')
            # print(for_collate)
            if for_collate:
                tbp_result = self._open_json(path = self.isolates[isolate]['collate'])
                raw_result = self._open_json(path = self.isolates[isolate]['profile'])
                _dict = self._infer_drugs(tbp_result = tbp_result,seq_id=isolate)
                # print(_dict)
                _dict = self._infer_dr_profile(res = _dict)
                _dict['Species'] = self._species(res = tbp_result, seq_id=isolate)
                _dict['Phylogenetic lineage'] = self._lineage(res = tbp_result, seq_id=isolate)
                _dict['Database version'] = self._db_version(res = raw_result)
                _dict['Median genome coverage'] = self._get_qc_feature(seq_id = isolate, res =tbp_result, val = 'median_coverage')
                _dict['Percentage reads mapped'] = self._get_qc_feature(seq_id = isolate,res = tbp_result, val = 'pct_reads_mapped')
                self._save_data(_data = [_dict], prefix = f"{isolate}/tbtamr")
                # self._save_csv(_)
                results.append(_dict)
            # else:
            #     logger.info(f'Collated results already exist for {isolate}.')
            # break
        
        logger.info(f"Saving collated data.")
        self._save_data(_data = results, prefix = "tbtamr")
        



#  headers for output Seq_ID, 
# Identification (WGS), Phylogenetic lineage,Predicted drug resist. summary:, Rifampicin,Isoniazid,Pyrazinamide,Ethambutol,Moxifloxacin,Amikacin,Cycloserine,Ethionamide,Para-aminosalicylic acid,Clofazimine,Delamanid,Bedaquiline,Linezolid,Database