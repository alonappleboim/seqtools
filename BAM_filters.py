from abc import ABCMeta, abstractmethod


class BAMFilter(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_name(self):
        """
        :return: return the name used by the user to invoke this filter
        """
        pass

    @abstractmethod
    def get_description(self):
        """
        :return: return a description of this filter
        """

    @abstractmethod
    def get_args(self):
        """
        :return: return a list of expected arguments: [(name, datatype, description),...]
        """
        pass

    @abstractmethod
    def filter(self, bamin, bamout, args, negate):
        """
        filter the data in (path) bamin into bamout. Sorting and indexing of output is
        performed by wrapping calls, so do not waste time on this.

        :param bamin:  path to input - an indexed and sorted bam file
        :param bamout: path to output
        :param args:   a dictionary of argname: value
        :param negate: boolean, whether to reverse the filter (i.e. filter complement)
        :return: False if an error occured, True otherwise
        """
        pass

class DuplicateFilter(BAMFilter):

    def get_name(self):
        return 'dup'

    def get_description(self):
        return 'remove artificially amplified reads'

    def get_description(self):
        return 'remove artificially amplified reads'

    def get_args(self):
        return [('type', str, 'the method of choice for removal: "umi&start", "umi&end", "umi&start&end"')]


class AlignmentQualityFilter(BAMFilter):

    def get_name(self):
        return 'qual'

    def get_description(self):
        return 'remove reads with low alignment quality (see bowtie/SAM documentation)'

    def get_args(self):
        return [('q', int, 'the method of choice for removal: "umi&start", "umi&end", "umi&start&end"')]
