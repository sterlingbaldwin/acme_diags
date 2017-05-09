import json, itertools, os, pdb
def create(json_file, custom_dir):
    """This function creates a json file for each variable, season, region
    and plev in the file specified json file."""

    f=open(json_file).read()
    plotsets = json.loads(f)

    keys_to_split = ["variables", "season", "region", "plevs"]
    output_file_names = []
    for plotset_id, plotset in plotsets.items():
        for plotset_case in plotset:
            basefile = plotset_case["case_id"]
            lists_to_split = []
            keys = []
            for key in keys_to_split:
                #pdb.set_trace()
                if key in plotset_case.keys():
                    if len(plotset_case[key]) is not 0:
                        keys += [key]
                        if type(plotset_case[key][0]) is float:
                            lists_to_split += [map(lambda n: '{}'.format(n), plotset_case[key] )]
                        else:
                            lists_to_split += [plotset_case[key]]

            for L in itertools.product(*lists_to_split):
                plotset_copy = plotset_case.copy()
                outfile = custom_dir + basefile + '_' + '_'.join(L) + '.json'
                for key, value in zip(keys, L):
                    plotset_copy[key] = [value]

                f = open( outfile, 'w' )
                json.dump({plotset_id: [plotset_copy] }, f)
                f.close()
                output_file_names += [outfile]
    return output_file_names
if __name__ == '__main__':
    outdir = './jsons/'
    outfiles = create('custom.json', outdir)

    cmdfile = "ps5_commands.sh"
    out = open(cmdfile, "w")

    #logger.info("List of commands is in: %s", cmdfile)
    #module use / usr / common / contrib / acme / modulefiles
    #module load uvcdat / batch

    dryrun = False
    sbatch_number = 1
    if dryrun:
        if sbatch_number > 0:
            print >> out, "#!/bin/bash"
            print >> out, """#SBATCH -p debug
#SBATCH -N %i
#SBATCH -D $HOME/acme_slurm/logs/
#SBATCH -t 00:30:00
#SBATCH -J plotset
#SBATCH -o set5_driver.o%%j

source activate 2.8
""" % (sbatch_number)

    for outfile in outfiles:
        print >> out, \
                 os.environ['HOME'] + '/acme-diags/acme_diags/plotset5/set5_driver.py ' + \
                 '-p ' + os.environ['HOME'] + '/acme_slurm/parameter.py ' + \
                 '-d '+ outfile
