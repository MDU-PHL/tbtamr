import pathlib, subprocess, os, psutil
from tbtamr.CustomLog import logger


class Tbtamr(object):
    """
    A base class for setting up tbtamr return a valid input object for subsequent steps
    """
    def __init__(self):
        

        
        self.one,self.five,self.fifteen = psutil.getloadavg()
        self.total_cores = os.cpu_count()
        self._cwd = pathlib.Path.cwd()
        # self.jobs = args.jobs # number of tbprofilers to run at a time
        # self.read1 = args.read1
        # self.read2 = args.read2
        # self.prefix = args.prefix
        # self.datafile = args.data
        # self.database = args.database

    # def _get_logger(self):
    #     logger =logging.getLogger(__name__) 
    #     logger.setLevel(logging.DEBUG)
    #     ch = logging.StreamHandler()
    #     ch.setLevel(logging.DEBUG)
    #     ch.setFormatter(CustomFormatter())
    #     fh = logging.FileHandler('tbtamr.log')
    #     fh.setLevel(logging.DEBUG)
    #     formatter = logging.Formatter('[%(levelname)s:%(asctime)s] %(message)s', datefmt='%Y%m-%d %I:%M:%S %p') 
    #     fh.setFormatter(formatter)
    #     logger.addHandler(ch) 
    #     logger.addHandler(fh)

    #     return logger

    def _run_cmd(self, cmd):
        
        """
        Use subprocess to run the command for tb-profiler
        """
        logger.info(f"Now running : {cmd}")
        p = subprocess.run(cmd, shell = True, capture_output = True, encoding = "utf-8")
        if p.returncode == 0:
            logger.info(f"{cmd} completed successfully. Will now move on to phenotype inferrence.")
            return True
        else:
            logger.critical(f"There appears to have been a problem with running {cmd}. The following error has been reported : \n {p.stderr}")
            raise SystemExit
    
    def _check_output_file(self, seq_id, step):

        wldcrd = f"{seq_id}/results/{seq_id}.results.json" if step == 'profile' else f"{seq_id}/tb-profiler_report.json"
        
        p = sorted(self._cwd.glob(wldcrd))
        if p != []:
            logger.info(f"{p[0]} has been found")
            return f"{p[0]}"
        else:
            return False

    def _check_output(self, isolates, step = 'profile'):
        
        for iso in isolates:
            present = self._check_output_file(seq_id= iso, step = step)
            if present:
                isolates[iso][step] = self._check_output_file(seq_id= iso, step = step)
            else:
                logger.critical(f"There seems to be a serious problem - files {iso} were not created. Please check logs and try again.")
                raise SystemExit
        logger.info(f"All files for step : {step} have been created.")
        return isolates

    def _set_threads(self, jobs):
        
        jobs = int(jobs)/2
        
        max_tbjob  = self.total_cores - max(self.one,self.five,self.fifteen)
        logger.info(f"The available cores is : {max_tbjob}")
        
        
        if int(jobs) == 0:
            logger.info(f"Number of TB-profiler jobs to run {max_tbjob}")
            return int(max_tbjob/2)
        elif int(jobs) <  max_tbjob/2:
            logger.info(f"Number of TB-profiler jobs to run {jobs}")
            return int(jobs)
        else:
            logger.info(f"Number of TB-profiler jobs to run {max_tbjob}")
            return int(max_tbjob/2)
            
    def _clean_cmd(self, path):

        cmd = f"rm -rf {path}"
        return cmd

    def _file_present(self, name):
        """
        check file is present
        :name is the path of the file to check
        """
        
        if name == "":
            return False
        elif pathlib.Path(name).exists():
            logger.info(f"Checking if file {name} exists")
            return True
        else:
            return False
