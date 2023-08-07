
import pathlib, subprocess, subprocess,re
from collections import namedtuple
from tbtamr.CustomLog import logger
from tbtamr.TbTamr import Tbtamr
from tbtamr.TbTamr_Utils import check

class RunProfiler(Tbtamr):
    """
    A class to run tbprofiler
    """
    def __init__(self, args):
        
        super().__init__()
        self.database = args.db
        self.input_data = args.input_data
        self.jobs = args.jobs
        self.keep = args.keep
        self.keep_bam = args.keep_bam
        self.exclude_not_reportable = args.exclude_not_reportable
        self.min_depth = args.min_depth
        self.prop_mtb = args.prop_mtb
        self.min_cov = args.min_cov
        
    def _get_isolates(self):

        isolates = {}
        if isinstance(self.input_data, dict):
            isolates[self.input_data['Seq_ID']] = {}
        else:
            try:
                with open(self.input_data,'r') as i:

                    lines = i.read().strip().split('\n')
                    for line in lines:
                        iso = line.split('\t')[0]
                        if iso not in isolates and iso != '':
                            isolates[iso] = {}
            except Exception as err:
                logger.critical(f"Something has gone wrong with your inputs. The following error was reported : {err}")
        return isolates
    
    def _remove(self, keep_bam = False, keep = False):

        if not keep:
            folders = ['vcf','bam'] if not keep_bam else ['vcf']
            for f in folders:
                p = pathlib.Path().resolve()
                cmd = self._clean_cmd(path = f"{p}/*/{f}")
                self._run_cmd(cmd=cmd)
        else:
            logger.info(f"Keeping all accessory data folders.")

    def _tidy_tbp(self):

        p = pathlib.Path().resolve()

        logger.info(f"Now tidying up")
        cmd = self._clean_cmd(path = f"{p}/*vcf*")

        self._run_cmd(cmd=cmd)

    def _single_cmd(self, input_data):
        cmd = f"tb-profiler profile --read1 {input_data['R1']} --read2 {input_data['R2']} {self.database} --prefix {input_data['Seq_ID']} --dir {input_data['Seq_ID']} --min_depth {self.min_depth} --no_trim --call_whole_genome --threads {self.jobs} >> {input_data['Seq_ID']}/tbprofiler.log 2>&1"
        return cmd

    def _batch_cmd(self, input_data):
        cmd = f"parallel --colsep '\\t' -j {self.jobs} 'tb-profiler profile --read1 {{2}} --read2 {{3}} {self.database} --prefix {{1}} --dir {{1}} --min_depth {self.min_depth} --no_trim --call_whole_genome --threads 1 >> {{1}}/tbprofiler.log 2>&1' :::: {input_data}"
        return cmd

    def _single_collate(self, input_data):
        cmd = f"tb-profiler collate -d {input_data['Seq_ID']}/results/ {self.database} -p {input_data['Seq_ID']}/tb-profiler_report --full --all_variants --mark_missing"
        return cmd

    def _batch_collate(self, input_data):
        cmd = f"parallel --colsep '\\t' -j {self.jobs} tb-profiler collate -d {{1}}/results/ {self.database} -p {{1}}/tb-profiler_report --full --all_variants --mark_missing :::: {input_data}"
        return cmd

    def _check_tbprofiler(self):
        
        return check()
        
        
    def _run(self):
        """
        run tbprofiler
        """
        
        if self._check_tbprofiler():
            
            isolates = self._get_isolates()
            if isolates != {}:
                logger.info(f"All check complete now running TB-profiler")
            else:
                logger.info(f"It seems that you have already run tbtamr on these sequences. Now exiting")
                raise SystemExit
        else:
            logger.critical(f"TB-profiler does not seem to be installed correctly. Please try again.")
            raise SystemExit
        # get list isolates
        
        # self._check_tbprofiler()
        cmd_profiler = self._single_cmd(input_data= self.input_data) if len(isolates) == 1 else self._batch_cmd(input_data= self.input_data)
        if self._run_cmd(cmd = cmd_profiler):
            isolates = self._check_output(isolates = isolates, step = 'profile')
            logger.info(f"Profiling was completed successfully, now collating results.")
            cmd_collate = self._single_collate(input_data= self.input_data) if len(isolates) == 1 else self._batch_collate(input_data=self.input_data)
            if self._run_cmd(cmd = cmd_collate):
                isolates = self._check_output(isolates = isolates, step = 'collate')
        # clean up
                self._tidy_tbp()
                self._remove(keep_bam=self.keep_bam, keep = self.keep)
                Input = namedtuple('Input', ['isolates',  'exclude_not_reportable','min_depth','min_cov', 'prop_mtb'])
                to_input = Input(isolates, self.exclude_not_reportable,self.min_depth, self.min_cov, self.prop_mtb) 
        
                return to_input
                
        
        logger.critical(f"Something seems to be wrong with your run of tbTAMR. Please try again.")


        
