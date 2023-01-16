import gzip


class OpenTextFile:
    """ Simple iterator to open standard text file formats
        Keyword Arguments
            text_file_path (str): Path to text file to iterate through
            separator (str): character used to split output string
            yield_indices (list, tuple): list or tuple of indices to yield if the entire line isn't needed
            yield_start (str): when str is observed lines are then yielded at the next line
            yield_end (str): when str is observed lines stop being yielded at the line
        Attributes:
            self.f (text.io): object to stream text/ byte lines
        """

    def __init__(self, text_file_path=None, separator='\t', yield_indices=None, yield_start=None, yield_end=None):
        assert isinstance(text_file_path, str)
        self.separator = separator
        self.yield_indices = yield_indices
        if yield_indices:
            assert isinstance(yield_indices, (list, tuple))
        self.yield_start = yield_start
        if self.yield_start:
            assert isinstance(yield_start, str)
        self.yield_end = yield_end
        if text_file_path.endswith(".gz"):
            self.file = gzip.open(text_file_path, 'rb')
        else:
            self.file = open(text_file_path, 'r')
        self.start_yield = True
        if self.yield_start:
            self.start_yield = False
        self.yield_stop = None

    def __iter__(self):
        # open file
        while True:
            # read line
            line = self.file.readline()
            # if line is blank break loop
            if not line or self.yield_stop:
                break
            # process line
            processed_line = self.process_line(line)
            # check if line should be yielded
            if self.check_yield(processed_line):
                yield processed_line

    def check_yield(self, line):
        # if yield identifier is None always return true
        if not self.start_yield:
            if self.yield_start in line[0]:
                # return false and on next iteration return True
                self.start_yield = True
            return False
        # if self.yield_end check for str in line and stop iteration at the line
        if self.yield_end:
            if self.yield_end in line[0]:
                self.yield_stop = True
                return False
        return True

    def process_line(self, line):
        # decode str if byte
        if isinstance(line, bytes):
            processed_line = line.decode('utf-8').strip().split(self.separator)
        else:
            processed_line = line.replace.strip().split(self.separator)
        # if yield_indices return processed line
        if self.yield_indices:
            if len(processed_line) >= len(self.yield_indices):
                return [processed_line[index] for index in self.yield_indices]
        return processed_line