import multiprocessing as mp
from multiprocessing import Pool


def process_list(list_to_process, process_function, n_threads=None, chunksize=5):
    """
    Process a function over elements in a list in parallel

    :param list_to_process: Will be deep-copied!
    :param process_function: Function to apply to each element. This should return the element after modification
    :param n_threads: Number of CPU threads, defaults to the number of cores
    :param chunksize: default 5
    :return: List of modified elements
    """
    if n_threads is None:
        n_threads = mp.cpu_count()

    pool = Pool(n_threads)
    result = pool.map(process_function, list_to_process, chunksize)
    pool.close()
    pool.join()

    return result


if __name__ == '__main__':
    pass
