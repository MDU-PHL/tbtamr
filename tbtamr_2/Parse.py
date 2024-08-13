import sys,gzip,pandas,pathlib,json,warnings,subprocess
from CustomLog import logger
warnings.simplefilter(action='ignore', category=FutureWarning)
pandas.set_option("mode.chained_assignment", None)

class Vcf(object):

    def __init__(self, 
                 vcf,
                 catalog,
                 catalog_config,
                 seq_id):
        self.vcf_file = vcf
        self.catalog = catalog
        self.config = self.get_config(pth = catalog_config)
        self.seq_id = seq_id

    def run_cmd(self, cmd) -> bool:

        logger.info(f"Now running {cmd}")

        proc = subprocess.run(cmd, shell = True, capture_output=True, encoding='utf-8')

        if proc.returncode == 0:
            logger.info(f"{proc.stdout}")
            return True
        else:
            logger.critical(f"{cmd} failed. The following error was encountered : {proc.stderr}")
            raise SystemExit
    def get_config(self, pth):

        if self.check_file(pth= pth):
            with open(pth, 'r') as j:

                return json.load(j)

    def check_file(self, pth) -> bool:

        if pth != "" and pathlib.Path(pth).exists():
            logger.info(f"{pth} exists.")
        else:
            logger.critical(f"{pth} does not exists or can't be accessed. Please check your inputs and try again.")
            raise SystemExit

        return True

    def check_gene(self, annot, genes) -> bool:
        g = annot.split('|')[3]
        if g in genes:
                return True        
        return False
    
    def calc_af(self, ao, dp):
        if ao != [] and dp != []:
            
            try:
                a = int(ao[0].split('=')[-1])
                d = int(dp[0].split('=')[-1])
                af = round(a/d, 2)
                return af
            except Exception as e:
                logger.critical(f"Something is wrong with the format of your vcf for {self.seq_id} the following error was encountered: {e} AO = {ao} and DP = {dp}")
                return 0
        return 0
    def check_lof(self,lof,gene,genes,vr, af):
        if not gene in genes:
            return vr,af
        else:
            ofint = [l for l in lof[0].split(',') if gene in l]
            if ofint != []:
                # print(ofint[0])
                a = ofint[0].split('|')[-1].strip(')')
                return f"{gene}_LoF",float(a)
            else:
                return vr,af
        

        
        

    def top_ann(self,nfo, genes ) -> dict:

        a = (i for i in nfo.split(';') if 'ANN' in i)
        ao = [i for i in nfo.split(';') if 'AO=' in i]
        dp = [i for i in nfo.split(';') if 'DP=' in i]
        lof = [i for i in nfo.split(';') if 'LOF' in i]
        af = self.calc_af(ao = ao , dp = dp)
        res = (r for r in a )
        annots = (v.split(',') for v in res)
        for annot in annots:
            # print(a)
            vals = (a.split('|') for a in annot if self.check_gene(annot = a, genes=genes))    
            # print(vals)
        cols = ['effect','gene','change_nuc','change_aa', 'variant', 'af']
        ix = [1,3,9,10]
        dt= []
        for record in vals:
            # print(record)
            rs = [record[r] for r in ix]
            vr = f"{rs[1]}_{rs[-1]}" if rs[-1] != '' else f"{rs[1]}_{rs[-2]}"
            if rs[0] in ['transcript_ablation','feature_ablation'] :
                vr = f"{rs[1]}_deletion"
            elif lof != []:
                vr,af = self.check_lof(lof = lof, gene = rs[1], genes = genes, vr = vr, af = af)
            rs.append(vr)
            rs.append(af)
            dt.append(dict(zip(cols, rs)))
            
        # print(dt)
        return dt

    def variant_generator(self,vcf_file, genes) -> list:

        all_lines = (line for line in gzip.open(vcf_file, 'r'))
        data = (s.decode().strip().split('\t') for s in all_lines if "##" not in s.decode() )
        cols = next(data)
        
        parsed_vcf = (dict(zip(cols, val)) for val in data)
        results = []
        for record in parsed_vcf:
            result = self.top_ann(nfo = record['INFO'], genes = genes)
            
            if result != []:
                results.extend(result)
        
        return results
    
    def create_output_dir(self,seq_id, force = False) -> bool:

        cmd = f"mkdir -p {seq_id}"
        logger.info(f"Will now create directory for {seq_id}")
        proc = self.run_cmd(cmd = cmd)
        if proc:
            return True
        
        return False
    
    def save_variants(self, df) -> bool:
        self.create_output_dir(seq_id=self.seq_id)
        pandas.DataFrame(df).to_csv(f'{self.seq_id}/{self.seq_id}_variants.csv', index = False)
        return True
    
    def get_catalog(self, catalog) -> pandas.DataFrame:
        return pandas.read_csv(catalog, dtype= str)
    
    def get_variant_data(self) -> list:
        if self.check_file(pth = self.vcf_file) and self.check_file(pth = self.catalog):
            try:
                catalog = self.get_catalog(catalog=self.catalog)
                genes = list(catalog[self.config['gene_name_col']].unique())
                variants = self.variant_generator(vcf_file = self.vcf_file, genes = genes)
                logger.info(f"Variants have been extracted from {self.vcf_file}")
                logger.info(f"Saving variants.")
                self.save_variants(df = variants)

                return variants
            except Exception as e:
                logger.critical(f"Something is wrong with your vcf format : {e}")
                raise SystemExit









