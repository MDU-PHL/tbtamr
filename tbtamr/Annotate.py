import sys,gzip,pandas,pathlib,json, subprocess, os,logging
from .CustomLog import logger
from .Utils import check_annotate


def check_file(pth) -> bool:

    if pth != "" and pathlib.Path(pth).exists():
        logger.info(f"{pth} exists.")
    else:
        logger.critical(f"{pth} does not exists or can't be accessed. Please check your inputs and try again.")
        raise SystemExit

    return True

def check_canannotate() -> bool:

    if check_annotate():
        return True

    else:
        logger.critical(f"Can't run annotation - you do not have snpEff installed properly. Exiting...")
        raise SystemExit

def run_cmd(cmd) -> bool:

        logger.info(f"Now running {cmd}")

        proc = subprocess.run(cmd, shell = True, capture_output=True, encoding='utf-8')

        if proc.returncode == 0:
            logger.info(f"{proc.stdout}")
            return True
        else:
            logger.critical(f"{cmd} failed. The following error was encountered : {proc.stderr} | {proc.stdout}")
            raise SystemExit

def wrangle_snpeffdb():

    paths = os.getenv('PATH').split(':')
    for pth in paths:
        d = pathlib.Path(pth.strip('bin'))
        # p = os.environ.get('CONDA_PREFIX').strip('bin')
        cfg = sorted(pathlib.Path(d,'share').glob('snpeff*/snpEff.config'))
        
        if cfg != []:
            return cfg[0]
    logger.critical(f"It seems that there is no config setup for snpEff. Please chack your installation and try again. Exiting....")
    raise SystemExit

def run_snpeff(vcf_file, seq_id) -> bool:

        spc = "Mycobacterium_tuberculosis_h37rv"
        cfg = wrangle_snpeffdb()
        nme = vcf_file.strip('.vcf.gz')
        snpeff =f"snpEff ann -dataDir . -c {cfg} -noLog -noStats {spc} {vcf_file} > {seq_id}/{seq_id}.annot.vcf"
        logger.info(f"Annotating vcf file")
        run_cmd(cmd=snpeff)
        run_cmd(cmd = f"bgzip {seq_id}/{seq_id}.annot.vcf")
        run_cmd(cmd = f"bcftools index {seq_id}/{seq_id}.annot.vcf.gz")

        return True

def create_output_dir(seq_id) -> bool:
  
    try:
        pathlib.Path(f"{seq_id}").mkdir(exist_ok=True)
        return True
    except:
        logger.critical(f"Something has gone wrong creating the folder for {seq_id}.")
        raise SystemExit
    

def annotate(vcf_file, seq_id):
    
    create_output_dir(seq_id=seq_id)
    fh = logging.FileHandler(f'{seq_id}/snpeff.log')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p') 
    fh.setFormatter(formatter)
    logger.addHandler(fh)

    if check_canannotate() and create_output_dir(seq_id = seq_id) and check_file(pth = vcf_file):
        run_snpeff(vcf_file= vcf_file, seq_id = seq_id)
        
        return f"{seq_id}/{seq_id}.annot.vcf.gz"
    else:
        logger.critical(f"Annotation cannot be performed - snpEff is not installed. Please check your installation and try again.")