import os,pathlib



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

    try:
        from pathogenprofiler import barcode, Vcf
        return True
    except:
        return False
