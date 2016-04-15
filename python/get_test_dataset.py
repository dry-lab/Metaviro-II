# encoding utf-8


def main(argv=None):
    logger.info("PREPROCESS")
    logger.info("Getting taxids...")

    with open("taxids.pickle", 'rb') as input:
        taxids = pickle.load(input)

    logger.info("TRAINING")

    random.seed()
    train_gis = [[]]*4

    cpt = 0

    output_files = [open("test_metaviro/test_data/test_data_%s" % pattern, 'w') for pattern in patterns]

    while cpt < 15000:
        for k in categories.keys():
            train_gis[k] = random_sampling(taxids[k])
                
        min_size = min([len(l) for l in train_gis])
        cpt += min_size
        
        for i in range(min_size):
            for cat in categories.keys():
                (gi, pos, read_size) = train_gis[cat][i]
                seq = get_NCBI_slice_sequence(gi, pos, pos+read_size-1)
                for pi in range(len(patterns)):
                    pattern = patterns[pi]
                    features = kmerize(seq, pattern, ref_read_size)
                    vw_feat = print_vw_features(int(gi), features, cat + 1)
                    output_files[pi].write(vw_feat)
        sleep(5)
        logger.info(cpt)

    for output_file in output_files:
        output_file.close()

    return 0

if __name__ == '__main__':
    sys.exit(main())


