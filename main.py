# expects two command line arguments:
# (1) random seed (use -1 for original seeds)
# (2) output directory

from distutils.dir_util import copy_tree
from shutil import copy
from sys import argv

import train.run_me
import analysis.prepare_data


if __name__ == "__main__":
    if len(argv) == 1 or argv[1] == "-1":
        random_seed = 234
        pipeline_seed = 20080808
    else:
        random_seed = int(argv[1])
        pipeline_seed = int(argv[1])

    print "########## Using seeds", random_seed, "and", pipeline_seed

    train.run_me.main(random_seed=random_seed, pipeline_seed=pipeline_seed)
    analysis.prepare_data.run()

    target_dir = "data/{}/{}_{}".format(argv[2], random_seed, pipeline_seed)

    print "########## Copying results to", target_dir

    copy_tree("../analysis/extracted/", target_dir)
    copy_tree("../_logs/p1000/pnet", target_dir)
