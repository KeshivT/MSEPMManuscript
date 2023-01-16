from TextIterator import OpenTextFile


class SeriesMatrixMetaExtractor:

    def __init__(self, series_matrix_path=None, sample_id=None, platform_id_filter='GPL13534',
                 phenotype_ids=None, stop_descriptor='!series_matrix_table_begin'):
        self.series_matrix = OpenTextFile(text_file_path=series_matrix_path, yield_end=stop_descriptor)
        self.sample_id = sample_id
        self.phenotype_ids = phenotype_ids
        self.phenotype_matrix = {}
        self.platform_filter = platform_id_filter
        self.series_description = {}
        self.series_matrix_raw_metadata = {}

    def get_meta_data(self):
        self.get_matrix_meta()
        self.experiment_id = self.series_matrix_raw_metadata['!Series_geo_accession'][0].strip('"')
        self.platform_id = self.series_matrix_raw_metadata['!Series_platform_id'][0].strip('"')
        if self.platform_id == self.platform_filter:
            self.data_processing = self.get_descriptive_lines('!Sample_data_processing')
            self.tissue_source = self.get_descriptive_lines('!Sample_source_name_ch1')
            self.sample_extraction = self.get_descriptive_lines('!Sample_extract_protocol_ch1')
            self.sample_ids = self.get_sample_ids
            self.get_matrix_meta()
            self.parse_metadata()

    def get_descriptive_lines(self, line_id):
        return [sample.strip('"') for sample in self.series_matrix_raw_metadata[line_id]]

    @property
    def get_sample_ids(self):
        sample_ids = [sample.strip('"') for sample in self.series_matrix_raw_metadata[self.sample_id]]
        for count, sample_id in enumerate(sample_ids):
            self.phenotype_matrix[sample_id.strip('"')] = {'experiment_id': self.experiment_id,
                                                           'platform_id': self.platform_id,
                                                           'tissue_source': self.tissue_source[count],
                                                           'data_processing': self.data_processing[count],
                                                           'sample_extraction': self.sample_extraction[count]}
        return sample_ids

    def get_matrix_meta(self):
        for count, line in enumerate(self.series_matrix):
            if len(line) > 1:
                if line[0] in self.series_matrix_raw_metadata.keys():
                    self.series_matrix_raw_metadata[f'{line[0]}___{count}'] = line[1:]
                else:
                    self.series_matrix_raw_metadata[f'{line[0]}'] = line[1:]

    def parse_metadata(self):
        for key, value in self.series_matrix_raw_metadata.items():
            if key.split('___')[0] in self.phenotype_ids:
                for phenotype, sample_id in zip(value, self.sample_ids):
                    phenotype_split = phenotype.strip(' "').split(':')
                    if len(phenotype_split) > 1:
                        phenotype_label = phenotype_split[0]
                        self.phenotype_matrix[sample_id][phenotype_label] = phenotype_split[1].lstrip()