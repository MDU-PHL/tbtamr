import sys,gzip,pandas,pathlib,json,warnings,subprocess,logging
from .CustomLog import logger
from .Annotate import annotate
warnings.simplefilter(action='ignore', category=FutureWarning)
pandas.set_option("mode.chained_assignment", None)

class Vcf(object):

    def __init__(self, 
                 vcf,
                 catalog,
                 catalog_config,
                 seq_id,
                 force):
        self.vcf_file = vcf
        self.catalog = catalog
        self.config = self.get_config(pth = catalog_config)
        self.seq_id = seq_id
        self.force = force

        self.create_output_dir(seq_id=self.seq_id, force= self.force)
        fh = logging.FileHandler(f'{seq_id}/tbtamr.log')
        fh.setLevel(logging.DEBUG)
        formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
        fh.setFormatter(formatter)
        logger.addHandler(fh)

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

    def check_data(self, vcf_file, _type = 'gzip') -> bool:
        
        if _type == 'gzip':
            all_lines = (line for line in gzip.open(vcf_file, 'r'))
            header = (s.decode() for s in all_lines if "##" in s.decode() )
        else:
            all_lines = (line for line in open(vcf_file, 'r'))
            header = (s for s in all_lines if "##"  in s )

        ao = False
        dp = False
        snpeff = False

        for record in header:
            if 'INFO=<ID=AO' in record:
                ao = True
            elif 'INFO=<ID=DP' in record:
                dp = True
            elif 'SnpEffCmd' in record:
                snpeff = True
        
        return ao and dp, snpeff
        
    def try_annotate(self, vcf_file) -> str:

        vcf_file = annotate(seq_id= self.seq_id, vcf_file = vcf_file)
        return vcf_file

    def get_data(self, vcf_file):

        try:
            all_lines = (line for line in gzip.open(vcf_file, 'r'))
            data = (s.decode().strip().split('\t') for s in all_lines if "##" not in s.decode() )
            _type = 'gzip'
        except:
            all_lines = (line for line in open(vcf_file, 'r'))
            data = (s.strip().split('\t') for s in all_lines if "##" not in s )
            _type = 'unzipped'
        
        dpths,annot = self.check_data(vcf_file = vcf_file, _type = _type)
        if dpths and annot:
            return data
        elif dpths and not annot:
            logger.warning(f"It looks like your vcf file is not correctly annotated - please wait while annotation is attempted.")
            vcf_file = self.try_annotate(vcf_file = vcf_file)
            all_lines = (line for line in gzip.open(vcf_file, 'r'))
            return (s.decode().strip().split('\t') for s in all_lines if "##" not in s.decode() )
        else:
            logger.critical(f"Something is wrong with your vcf file format. Please read documentation and try again. Exiting..")
            raise SystemExit

    def variant_generator(self,vcf_file, genes) -> list:

        
        data = self.get_data(vcf_file = vcf_file)
        cols = next(data)
        
        parsed_vcf = (dict(zip(cols, val)) for val in data)
        results = []
        for record in parsed_vcf:
            result = self.top_ann(nfo = record['INFO'], genes = genes)
            
            if result != []:
                results.extend(result)
        
        return results
    
    def create_output_dir(self,seq_id, force) -> bool:

        logger.info(f"Will now create directory for {seq_id}")
        # ex = not force
        # pri
        msg = f"Something has gone wrong creating the folder for {seq_id}." if force else f"Folder for {seq_id} exists. Please use --force if you would like to override any previous results."
        try:
            pathlib.Path(f"{seq_id}").mkdir(exist_ok=force)
            return True
        except:
            logger.critical(msg)
            raise SystemExit
        
    
    def save_variants(self, df) -> bool:
        self.create_output_dir(seq_id=self.seq_id, force = self.force)
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
                logger.critical(f"Something is wrong with your vcf format : {e} If you are unsure you can run try to run tbtamr in annotate mode (providing you have snpEff installed) or fqtovcf mode (providing that mutamr is installed).")
                raise SystemExit









