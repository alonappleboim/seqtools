
class TTScount(Aggregator):
    name = "orf_counter"
    description = "count reads that fall in certain window around TSS/TTS"
    args = {'f': ('TTS', 'around which feature to count, TSS/TTS'),
            'w': ([-100,100], 'window around feature in which reads are summed'),
            't': ('cov', "which kind of aggregation: 5',3',cov")}
    suffixes = {'.w.bw': ('bigwig', "watson strand coverage"),
                '.c.bw': ('bigwig', "crick strand coverage")}
    out_folder = "counter"

    def analyze(self, files):
        tmp = self.tmp_dir + os.sep + os.path.split(files['.w.bw'])[1].split(os.extsep)[0] + '.cov.bed.tmp'
        def handle_strand(char, fout):
            bedcmd = "bedtools genomecov -ibam %s -g %s -bg -strand %s"
            bed = sp.Popen(sh.split(bedcmd % (files['bam_in'], SGLP, STRANDS[char])), stdout=sp.PIPE)
            sbed = sp.Popen(sh.split("sort -k1,1 -k2,2n"), stdin=bed.stdout, stdout=open(tmp, 'w'))
            sbed.wait()
            bw = sp.Popen([BG2W_EXEC, tmp, SGLP, fout])
            bw.wait()
            os.remove(tmp)
        handle_strand('w', files['.w.bw'])
        handle_strand('c', files['.c.bw'])
        return {}

    def aggregate(self, filemap):
        pass


