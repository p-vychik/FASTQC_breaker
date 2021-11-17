import math
import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
def basic_statistics(sequences, fastq_name):

    # words to replace in html report
    words = {
        'BASIC_STATISTICS_filename':        fastq_name,
        'BASIC_STATISTICS_file_type':       'set here filetype',
        'BASIC_STATISTICS_encoding':        'set here encoding',
        'BASIC_STATISTICS_total_sequences': 'set here total seqs num',
        'BASIC_STATISTICS_seq_poor_qual':   'set here poor qual num',
        'BASIC_STATISTICS_seq_len':         'set here seq len',
        'BASIC_STATISTICS_gc':              'set here %GC'
    }
    return words


# -----------------------------------------------------------------------------
# create image and return status
def per_base_sequence_quality(sequences, fastq_name, imgname):

    # Creating dataset
    # go by positions
    all_quals = []
    bp_cnt = 0  # number of boxplot
    pos = 0     # position in sequence
    x_step = 1
    x_labels = []
    mean_line_x = []
    mean_line_y = []

    # get transpose matrix ant go throuhg line by line
    qual_for_base_mat = sequences.get_transpose_qual_matrix()
    xlen = qual_for_base_mat.shape[0]
    while 1:
        if xlen > 60 and bp_cnt == 9:
            # make groups - only 50 boxplots on X axis for x > 10
            x_step = math.ceil((xlen - 10) / 50)
        # combine x_step lines together
        quals = np.zeros(0, dtype=np.int8)
        for k in range(x_step - 0):
            idx = pos + 0 + k
            if idx >= xlen:
                break
            arr = qual_for_base_mat[idx, :]
            # arr can have '-1' elements at the end, remove them
            arr = arr[arr >= 0]
            if idx < xlen:
                quals = np.concatenate((quals, arr))
        # save
        all_quals.append(quals)
        mean_line_x.append(bp_cnt + 1)
        mean_line_y.append(sum(quals) / len(quals))
        # add pretty x-labels, make low density to be readable
        if x_step == 1:
            if bp_cnt < 9:
                x_labels.append(str(pos + 1))
            else:
                if ((pos) % 2) == 0:
                    x_labels.append(str(pos + 1))
                else:
                    x_labels.append("")  # empty label
        else:  # xstep > 1
            bp_offset = 2
            density = 4                  # print label each 'density' boxplot
            if xlen > 100:
                density = 5
            if ((bp_cnt + bp_offset) % density) == 0:
                # pair label
                first_num = pos + 1
                second_num = min(pos + x_step, xlen)
                x_labels.append(str(first_num) + '-' + str(second_num))
            else:
                x_labels.append("")      # empty label
        pos += x_step
        bp_cnt += 1
        if pos >= xlen:
            break

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    # creating axes instance
    bp = ax.boxplot(all_quals, patch_artist=True, showfliers=False,
                    whis=[10, 90],
                    boxprops=dict(facecolor='yellow', alpha=0.7))

    # to define ERROR or WARNING observ min_quartile and min_mean
    whiskers = [item.get_ydata() for item in bp['whiskers']]
    lower_whiskers = whiskers[::2]
    lower_quartiles = [item[0] for item in lower_whiskers]
    min_quartile = min(lower_quartiles)
    min_mean = min(mean_line_y)
    status = "good"
    if min_quartile < 10 or min_mean < 25:
        status = 'warning'
    if min_quartile < 5 or min_mean < 20:
        status = 'fail'

    # removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # make grid
    max_y = max(sequences.max_qual(), 30)
    plt.grid(axis='x')
    plt.axhspan(0, 20, color='red', alpha=0.2, zorder=0)
    plt.axhspan(20, 28, color='yellow', alpha=0.2, zorder=0)
    plt.axhspan(28, max_y, color='green', alpha=0.2, zorder=0)
    ax.set_xticklabels(x_labels)
    ax.set_ylim([0, sequences.max_qual()])

    # add title
    ax.set(xlabel='Position in read (bp)', ylabel='Quality score',
           title='Quality scores across all bases')

    # add mean line
    ax.plot(mean_line_x, mean_line_y)

    fig.savefig(imgname)
    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_sequence_quality_scores(sequences, fastq_name, imgname):

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)

    # get dataset
    meanq_arr, min_mq, max_mq = sequences.get_meanq_arr()

    # removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # make grid
    plt.grid(axis='x')
    plt.grid(axis='y')
    plt.axhspan(min_mq, max_mq, color='black', alpha=0.01, zorder=0)
    ax.set_xticks(range(min_mq, max_mq, 1))
    ax.set_xlim([min_mq, max_mq])
    ax.set_ylim([0, max(meanq_arr)])

    # add title
    ax.set(xlabel='Mean Sequence Quality', ylabel='Number of reads',
           title='Quality score destribution oer all sequences ' +
                 '(average quality per read)')

    # add mean line
    ax.plot(meanq_arr, color='red')

    # define status
    status = 'good'
    most_freq_meanq = np.where(meanq_arr == np.amax(meanq_arr))[0]
    if most_freq_meanq < 27:
        status = 'warning'
    if most_freq_meanq < 20:
        status = 'fail'

    fig.savefig(imgname)
    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_base_sequence_content(sequences, fastq_name, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='per_base_sequence_content')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'

    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_sequence_gc_content(sequences, fastq_name, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='per_sequence_gc_content')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_base_n_content(sequences, fastq_name, imgname):
    
    # set 2D array of seq
    reads = sequences.seq_mat

    # create empty list for future work
    list_of_N_content = [0] * reads.shape[1]
    
    # check if element of read has N base
    for i in range(reads.shape[0]):
        for j in range(reads.shape[1]):
            if reads[i][j] == ord('N'):
                list_of_N_content[j] += 1


    read_len = reads.shape[1]
    
    # translate to %
    list_of_N_content = [((x / read_len) * 100) for x in list_of_N_content]

    # create plot for n base
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.plot(list_of_N_content, color = 'red')
    ax.set(title='N content across all bases', xlabel = 'Position in read(bp)')
    ax.set_ylim(-2.0,100.0)
    ax.text(read_len-3, 97, r'%N', color = 'red')
    ax.grid(color = 'black', linestyle = '--', linewidth = 0.5)
    fig.savefig(imgname)

    # define report status
    status = 'good'

    if max(list_of_N_content) > 5:
        status = 'warning'
    if max(list_of_N_content) > 20:
        status = 'fail'
     
    return status


# -----------------------------------------------------------------------------
# create image and return status
def sequence_length_distribution(sequences, fastq_name, imgname):

     # set 2D array of seq
    reads = sequences.seq_mat

    # create empty list for future work
    des_of_length = np.arange(1, reads.shape[1]+1)
    length_count = np.zeros(reads.shape[1])

    for i in range(reads.shape[0]):
        
        line_read = reads[i,:]
        # remove empty elements
        line_read = line_read[line_read >= 0]
        # add plus 1 in position of length
        length_count[line_read.size-1] += 1

    where_zero = np.where(length_count == 0.0)
    des_of_length = np.delete(des_of_length,where_zero)
    length_count = np.delete(length_count,where_zero) 

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    plt.plot(des_of_length,length_count, color = 'red')
    plt.grid(axis = 'y')
    plt.title('Distribution of sequence lengths over all sequences')
    plt.xlabel('Sequence Length(bp)')
    plt.text(max(des_of_length)-3,max(length_count)+10,'Sequence Length', color = 'red')


    plt.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def sequence_duplication_levels(sequences, fastq_name, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='sequence_duplication_levels')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def overrepresented_sequences(sequences, fastq_name, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='overrepresented_sequences')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def adapter_content(sequences, fastq_name, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='adapter_content')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status
