from importlib.resources import path
import pathlib, pandas,re, datetime, subprocess, os, logging,subprocess,collections, psutil
from tbtamr.CustomLog import logger
from tbtamr.TbTamr import Tbtamr


class AmrSetup(Tbtamr):

    def __init__(self, args):
        super().__init__()
        self.jobs = args.jobs # number of tbprofilers to run at a time
        self.read1 = args.read1
        self.read2 = args.read2
        self.prefix = args.prefix
        self.datafile = args.input
        self.database = args.database
        self.keep = args.keep
        self.keep_bam = args.keep_bam
        # self.logger = self._get_logger()

    def _check_prefix(self):
        """
        If running type is not batch, then check that prefix is present and a string
        """

        if self.prefix == '':
            logger.critical(f"You must supply a sample or sequence id.")
            raise SystemExit
        else:
            return True

    def _check_reads(self, read):

        if not self._file_present(read):
            logger.critical(f"{read} is not a valid file path. Please check your input and try again.")
            raise SystemExit
        return True

    def _link_dirs(self, r1, r2, p):

        r1_source = pathlib.Path(r1).resolve()
        r2_source = pathlib.Path(r2).resolve()
        _dir = pathlib.Path(p).resolve()

        if not _dir.exists():
            _dir.mkdir(parents=True, exist_ok=True)
        
        
        r1_target = _dir / f"R1.fq.gz"
        r2_target = _dir / f"R2.fq.gz"

        r1_target.symlink_to(r1_source) if not r1_target.exists() else logger.info(f"{r1_target} exists - nothing else to be done.")
        r2_target.symlink_to(r2_source) if not r2_target.exists() else logger.info(f"{r2_target} exists - nothing else to be done.")

        return f"{p}\t{r1_target}\t{r2_target}"

    def _write_inputfile(self, input_lines):

        pth = pathlib.Path('input.tab').resolve()
        logger.info(f"Now writing input file.")
        pth.write_text('\n'.join(input_lines))

        return f"{pth}"

    def _input_files(self, jobs):
        """
        Ensure that the files (either contigs or amrfinder output) exist and return running type
        """
        
        input_lines = []
        if self.datafile != '':
            logger.info(f"Checking that the input data is present.")
            
            try:

                with open(self.datafile, 'r') as c:
                    data = c.read().strip().split('\n')
                    for line in data:
                        row = line.split('\t')
                        iso = self._check_output_file(seq_id=row[0], step = 'initial')
                        if iso == False:
                            if self._check_reads(row[1]) and self._check_reads(row[2]):
                                logger.info(f"Reads for {row[0]} have been found")
                                logger.info(f"Will now setup sample directory for {row[0]}.")
                                input_line = self._link_dirs(r1 = row[1], r2 = row[2], p =row[0] )
                                input_lines.append(input_line)
                        else:
                            logger.info(f"It seems that results already exist for {iso}.")
            except Exception as e:
                logger.critical(f"Something has gone wrong with your input file. Please try again {e}.")
                raise SystemExit
        elif self.datafile == '' and self._check_reads(self.read1) and self._check_reads(self.read2):
            self._check_prefix()
            logger.info(f"Will now setup sample directory for {self.prefix}.")
            input_line = self._link_dirs(r1 = self.read1, r2 = self.read2, p =self.prefix)
            input_lines.append(input_line)
            # jobs = self._set_threads(jobs=jobs)
            logger.info(f"Reads for {self.prefix} are present. tbTAMR can proceed.")
        else:
            logger.critical(f"Something has gone wrong with your inputs. Please try again.")
            raise SystemExit

        self.input_file = self._write_inputfile(input_lines = input_lines)
        jobs = self._set_threads(jobs=jobs)
        return jobs
   

    def _setup(self):
        
        # check that inputs are correct and files are present
        jobs = self._input_files(jobs = self.jobs)
        # check that prefix is present (if needed)
        
        if self.keep and not self.keep_bam:
            logger.info(f"Keep all accessory files.")
        elif self.keep_bam and not self.keep:
            logger.info(f"Keeping only bam files")
        elif self.keep and self.keep_bam:
            logger.info(f"You have elected to keep all and keep bam - we assume you mean you want to keep everything.")
        elif not self.keep and not self.keep_bam:
            logger.info(f"You have not decided to keep accesory files. All accessory and intermediate files will be removed following successful completion of tbTAMR.")
        
        Data = collections.namedtuple('Data', ['input_file', 'jobs', 'db', 'keep','keep_bam'])
        input_data = Data(self.input_file, jobs, self.database, self.keep, self.keep_bam)
        
        return input_data
    

