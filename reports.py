import numpy as np
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------------
def make_report_basic_statistics(sequences, filename):

    # words to replace in html report
    words = {
        'BASIC_STATISTICS_filename':        filename,
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
def make_report_per_base_sequence_quality(sequences, imgname):

    # Creating dataset
    # go by positions
    all_quals = []
    cnt = 0  # number of data set for one boxplot
    pos = 0  # position in sequence
    x_step = 1
    x_labels = []
    mean_line_x = []
    mean_line_y = []

    # get transpose matrix ant go throuhg line by line
    qual_for_base_mat = sequences.get_transpose_qual_matrix()
    xlen = qual_for_base_mat.shape[0]
    while 1:
        if xlen > 60 and cnt == 9:
            x_step = 2
        # combine x_step lines together
        quals = np.zeros(0, dtype=np.int8)
        for k in range(x_step - 0):
            idx = pos + 0 + k
            arr = qual_for_base_mat[idx, :]
            # arr can have '-1' elements at the end, remove them
            arr = arr[arr >= 0]
            if idx < xlen:
                quals = np.concatenate((quals, arr))
        # save
        all_quals.append(quals)
        mean_line_x.append(cnt + 1)
        mean_line_y.append(sum(quals) / len(quals))
        # add pretty x-labels, make low density to be readable
        if x_step == 1:
            if cnt < 9:
                x_labels.append(str(pos + 1))
            else:
                if ((pos) % 2) == 0:
                    x_labels.append(str(pos + 1))
                else:
                    x_labels.append("")  # empty label
        else:  # xstep > 1
            if ((pos+1) % 6) == 0:
                # pair label
                x_labels.append(str(pos + 1) + '-' + str(pos + 2))
            else:
                x_labels.append("")  # empty label
        pos += x_step
        cnt += 1
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
def make_report_per_sequence_quality_scores(sequences, imgname):

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
def make_report_per_base_sequence_content(sequences, imgname):
    lines = sequences.seq_mat
    a_count = np.count_nonzero(lines == ord('A'), axis=0)
    g_count = np.count_nonzero(lines == ord('G'), axis=0)
    c_count = np.count_nonzero(lines == ord('C'), axis=0)
    t_count = np.count_nonzero(lines == ord('T'), axis=0)
    percents = np.vstack([g_count*100/lines.shape[0], a_count*100/lines.shape[0], t_count*100/lines.shape[0], c_count*100/lines.shape[0]])

    # image creation
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='Sequence content across all bases', xlabel='Position in read (bp)')
    xs = range(sequences.max_len)  # max seq len in class sequences
    ys1 = percents[0]  # g
    ys2 = percents[1]  # a
    ys3 = percents[2]  # t
    ys4 = percents[3]  # c

    # add pretty x-labels, make low density to be readable
    cnt = 0
    pos = 0  # position in sequence
    x_step = 1
    x_labels = []  # list of labels
    xlen = lines.shape[1]  # length of reads
    while pos < xlen:  # iter by position
        if xlen > 60 and cnt == 9:
            x_step = 2
        if x_step == 1:
            if cnt < 9:
                x_labels.append(str(pos + 1))
            else:
                if ((pos) % 2) == 0:
                    x_labels.append(str(pos + 1))
                else:
                    x_labels.append("")  # empty label
        else:  # xstep > 1
            if ((pos + 1) % 6) == 0:
            # pair label
                x_labels.append(str(pos + 1) + '-' + str(pos + 2))
            else:
                x_labels.append("")  # empty label
        pos += 1
        cnt += 1

    ax.plot(xs, ys1, label='%G', color='r')
    ax.plot(xs, ys2, label='%A', color='blue')
    ax.plot(xs, ys3, label='%T', color='black')
    ax.plot(xs, ys4, label='%C', color='lime')
    ax.set_xlim([1, lines.shape[1] + 1])
    ax.set_xticks(range(1, lines.shape[1]+1, 1))
    ax.set_xticklabels(x_labels)
    ax.set_yticks(range(0, 110, 10))
    ax.set_ylim([0.0, 100.0])
    plt.legend()
    ax.grid()

    fig.savefig(imgname)

    # define report status
    GC_diff = percents[0] - percents[3]
    AT_diff = percents[1] - percents[2]
    status = 'good'
    for i in range((GC_diff.shape[0])):
        if 20 > np.sqrt(GC_diff[i]**2) >= 10 or 20 > np.sqrt(AT_diff[i]**2) >= 10:
            status = 'warning'
        elif np.sqrt(GC_diff[i]**2) >= 20 or np.sqrt(AT_diff[i]**2) >= 20:
            status = 'fail'
    print(status)
    return status


# -----------------------------------------------------------------------------
# create image and return status
def make_report_per_sequence_gc_content(sequences, imgname):
    # preparing data for graphics
    lines = sequences.seq_mat
    g_count = np.count_nonzero(lines == ord('G'), axis=1)
    c_count = np.count_nonzero(lines == ord('C'), axis=1)
    gc_count = g_count * 100 / lines.shape[1] + c_count * 100 / lines.shape[1]

    gc_mean = np.mean(gc_count)
    gc_std = np.std(gc_count)

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
def make_report_per_base_n_content(sequences, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='per_base_n_content(')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def make_report_sequence_length_distribution(sequences, imgname):

    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='sequence_length_distribution')
    fig.savefig(imgname)

    # define report status
    status = 'good'
    status = 'warning'
    status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def make_report_sequence_duplication_levels(sequences, imgname):

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
def make_report_overrepresented_sequences(sequences, imgname):

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
def make_report_adapter_content(sequences, imgname):

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