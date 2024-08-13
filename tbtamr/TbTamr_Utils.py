
import pathlib, subprocess, re,json
from collections import namedtuple
from tbtamr.CustomLog import logger
# from tbtamr.RunProfiler import RunProfiler

# version_pat_3 = 

def _run_cmd(cmd):
    p = subprocess.run(f"{cmd}", capture_output=True, encoding = "utf-8", shell = True)
    return p  

def version_pattern():
    return re.compile(r'\bv?(?P<major>[0-9]+)\.(?P<minor>[0-9]+)(?:\.(?P<release>[0-9]+)*)?(?:\.(?P<build>[0-9]+)*)?\b')

def _run_check(_input, sft):
    p = _run_cmd(cmd=_input)
    _str = p.stdout
    v = version_pattern().search(_str.strip())
    if v:
        v = v.group(0)
        logger.info(f"{sft} version {v} detected.")              
    elif sft == "snpEff" and p.returncode == 0:
        logger.info(f"{' '.join(p.stdout.split())} detected.") 
    elif sft == "tb-profiler":
        logger.warning(f"It seems that {sft} is not installed correctly.")
        res = input("Would you like to install it now? (y/n): ")
        if res.lower() == 'y':
            _install_tbprofiler()
        elif res.lower() == 'n':
            logger.critical(f"{sft} is required - please run 'tbtamr setup' now to proceed.")
            raise SystemExit
    else:
        logger.critical(f"{sft} is not installed.")
        raise SystemExit

def check_deps():

    with open(f"{pathlib.Path(__file__).parent / 'dep_config.json'}", "r") as j:

        software = json.load(j)
        logger.info(f"Checking that dependencies are installed and recording version.")
        
        for sft in software:
            _run_check(_input =software[sft], sft = sft)
        
        return True
   
def _install_tbprofiler():

    pp = _run_cmd(cmd="pip3 install git+https://github.com/MDU-PHL/pathogen-profiler")
    if pp.returncode == 0:
        logger.info(f"pathogen-profiler from https://github.com/MDU-PHL/ installed.")
        tbp = _run_cmd(cmd = "pip3 install git+https://github.com/MDU-PHL/TBProfiler")
        if tbp.returncode == 0:
            logger.info(f"TB-profiler https://github.com/MDU-PHL/ installed. You should now be good to go.")
            

def check():

    check_deps()
    return True
    
def install():
    _install_tbprofiler()