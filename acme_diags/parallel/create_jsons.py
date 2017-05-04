import json, itertools, pdb
def create(json_file, custom_dir):
    """This function creates a json file for each variable ans season
    in the file specified"""

    f=open(json_file).read()
    plotsets = json.loads(f)

    keys_to_split = ["variables", "season", "region", "plevs"]
    template = []
    for plotset_id, plotset in plotsets.items():
        for plotset_case in plotset:
            template += [plotset_case]
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
                outfile = basefile + '_' + '_'.join(L) + '.json'
                for key, value in zip(keys, L):
                    plotset_copy[key] = [value]

                f = open(custom_dir + outfile, 'w')
                json.dump({plotset_id: [plotset_copy] }, f)
                f.close()

if __name__ == '__main__':
    create('custom.json', './slurm_jsons/')

    if opts["dryrun"]:
        fnm = os.path.join(outpath, "metadiags_commands.sh")
        dryrun = open(fnm, "w")

        logger.info("List of commands is in: %s",fnm)

        if opts["sbatch"] > 0:
            print >> dryrun, "#!/bin/bash"
            print >> dryrun, """#SBATCH -p debug
#SBATCH -N %i
#SBATCH -t 00:30:00
#SBATCH -J plotset
#SBATCH -o set5_driver.o%%j

module use /usr/common/contrib/acme/modulefiles
module load uvcdat/batch
""" % (opts["sbatch"])
    else:
        dryrun = False
