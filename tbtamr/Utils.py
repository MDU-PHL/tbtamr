import os,pathlib
from .CustomLog import logger


def check_install(tool) -> bool:

    paths = os.getenv('PATH').split(':')
    for pth in paths:
        d = pathlib.Path(pth)
        tl = d /f"{tool}"
        if tl.exists():
            return True
    return False

def check_annotate():

    if check_install(tool = 'snpEff'):
        return True
    else:
        return False

def check_mutamr():

    if check_install(tool = 'mutamr'):
        return True
    else:
        return False


def check_lineage():
    logger.info(f"Will check if lineage can be run.")
    try:
        from pathogenprofiler import barcode, Vcf
        logger.info(f"Lineage calling can be undertaken.")
        return True
    except:
        logger.warning(f"Lineage calling cannot be undertaken - pathogenprofiler needs to be installed.")
        return False
