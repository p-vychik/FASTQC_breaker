import math
import numpy as np
import matplotlib.pyplot as plt
import hashlib
import plotly.graph_objects as go
import pandas as pd
from scipy import stats



# -----------------------------------------------------------------------------
def basic_statistics(sequences, fastq_name):
    # words to replace in html report
    words = {
        'BASIC_STATISTICS_filename': fastq_name,
        'BASIC_STATISTICS_file_type': 'set here filetype',
        'BASIC_STATISTICS_encoding': 'set here encoding',
        'BASIC_STATISTICS_total_sequences': 'set here total seqs num',
        'BASIC_STATISTICS_seq_poor_qual': 'set here poor qual num',
        'BASIC_STATISTICS_seq_len': 'set here seq len',
        'BASIC_STATISTICS_gc': 'set here %GC'
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
    lines = sequences.seq_mat

    a_count = np.count_nonzero(lines == ord('A'), axis=0)
    g_count = np.count_nonzero(lines == ord('G'), axis=0)
    c_count = np.count_nonzero(lines == ord('C'), axis=0)
    t_count = np.count_nonzero(lines == ord('T'), axis=0)
    reads_length = np.count_nonzero(lines != -1, axis=0)
    percents = np.vstack([g_count * 100 / reads_length, a_count * 100 / reads_length,
                          t_count * 100 / reads_length, c_count * 100 / reads_length])
    # image creation
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='Sequence content across all bases', xlabel='Position in read (bp)')
    xs = range(sequences.max_len)  # max seq len in class sequences
    ys1 = percents[0]  # g
    ys2 = percents[1]  # a
    ys3 = percents[2]  # t
    ys4 = percents[3]  # c

    x_labels = []
    for i in range(lines.shape[1]):
        if (i < 101 and i % 5 == 0) or (i > 100 and i % 10 == 0):
            x_labels.append(i)
        else:
            x_labels.append('')

    ax.plot(xs, ys1, label='%G', color='r')
    ax.plot(xs, ys2, label='%A', color='blue')
    ax.plot(xs, ys3, label='%T', color='black')
    ax.plot(xs, ys4, label='%C', color='lime')
    ax.set_xlim([1, lines.shape[1] + 1])
    ax.set_xticks(range(1, lines.shape[1] + 1, 1))
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
        if 20 >= np.sqrt(GC_diff[i] ** 2) > 10 or 20 >= np.sqrt(AT_diff[i] ** 2) > 10:
            status = 'warning'
        elif np.sqrt(GC_diff[i] ** 2) > 20 or np.sqrt(AT_diff[i] ** 2) > 20:
            status = 'fail'

    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_base_gc_content(sequences, fastq_name, imgname):
    lines = sequences.seq_mat

    # preparing data
    g_count = np.count_nonzero(lines == ord('G'), axis=0)
    c_count = np.count_nonzero(lines == ord('C'), axis=0)
    reads_length = np.count_nonzero(lines != -1, axis=0)
    gc_count = (g_count + c_count) * 100 / reads_length

    # image creation
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='GC content across all bases')
    xs = range(sequences.max_len)  # max seq len in class sequences

    x_labels = []
    for i in range(lines.shape[1]):
        if (i < 101 and i % 5 == 0) or (i > 100 and i % 10 == 0):
            x_labels.append(i)
        else:
            x_labels.append('')

    ax.plot(xs, gc_count, label='%GC', color='r')
    ax.set_xlim([1, lines.shape[1] + 1])
    ax.set_xticks(range(1, lines.shape[1] + 1, 1))
    ax.set_xticklabels(x_labels)
    ax.set_yticks(range(0, 110, 10))
    ax.set_ylim([0.0, 100.0])
    plt.legend()
    ax.grid()

    fig.savefig(imgname)

    # define status
    gc_mean = np.mean(gc_count)
    status = 'good'
    for i in gc_count:
        if 5 <= np.abs(i - gc_mean) <= 10:
            status = 'warning'
        elif np.abs(i - gc_mean) < 30:
            status = 'fail'
        return status


# -----------------------------------------------------------------------------
# create image and return status
def per_sequence_gc_content(sequences, fastq_name, imgname):
    lines = sequences.seq_mat

    # preparing data
    g_count = np.count_nonzero(lines == ord('G'), axis=1)
    c_count = np.count_nonzero(lines == ord('C'), axis=1)
    reads_length = np.count_nonzero(lines != -1, axis=1)
    gc_count = (g_count + c_count) * 100 / reads_length

    gc_mean = np.median(gc_count)
    gc_std = np.std(gc_count)
    gc_theor = np.random.normal(loc=gc_mean, scale=gc_std, size=lines.shape[0])

    # image creation
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='GC distribution over all sequences')

    bins = np.linspace(0, lines.shape[1], lines.shape[1])
    gc_count_kernel = stats.gaussian_kde(gc_count)
    gc_count_curve = gc_count_kernel(bins) * gc_count.shape[0]
    gc_theor_kernel = stats.gaussian_kde(gc_theor)
    gc_theor_curve = gc_theor_kernel(bins) * gc_count.shape[0]

    ax.plot(bins, gc_count_curve, color="red", label='GC count per read')
    ax.plot(bins, gc_theor_curve, color="blue", label='Theoretical Distribution')
    plt.legend()

    ax.set_xticks(range(0, lines.shape[1], 1))
    x_labels = []
    for i in range(lines.shape[1]):
        if (i < 101 and i % 5 == 0) or (i > 100 and i % 10 == 0):
            x_labels.append(i)
        else:
            x_labels.append('')
    ax.set_xticklabels(x_labels)
    plt.xlabel("Mean GC content (%)")
    plt.grid()
    fig.savefig(imgname)

    # define report status
    dev_percent = 0
    theor = 0
    for i in range(gc_count.shape[0]):
        theor += gc_count[i]
        dev_percent += np.abs(gc_count[i] - gc_theor[i])
    dev_percent = (dev_percent / theor) * 100
    status = 'good'
    if 15 <= dev_percent < 30:
        status = 'warning'
    elif dev_percent >= 30:
        status = 'fail'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def per_base_n_content(sequences, fastq_name, imgname):
    # set 2D array of seq
    reads = sequences.seq_mat
    # create list with N bases
    list_of_N_content = np.count_nonzero(reads == ord('N'), axis=0)
    # translate to %
    list_of_N_content = list_of_N_content * 100 / reads.shape[1]
    # create plot for n base
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.plot(list_of_N_content, color='red')
    ax.set(title='N content across all bases', xlabel='Position in read(bp)')
    ax.set_ylim(-2.0, 100.0)
    ax.text(reads.shape[1] - 3, 97, r'%N', color='red')
    ax.grid(color='black', linestyle='--', linewidth=0.5)
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
    length_count = [0] * reads.shape[0]

    for i in range(reads.shape[0]):
        line_read = reads[i, :]
        # remove empty elements
        line_read = line_read[line_read >= 0]
        # add plus 1 in position of length
        length_count[i] += line_read.size

    mi_len = min(length_count)
    ma_len = max(length_count)
    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.hist(length_count, color='red', align='left')
    ax.set_xticks([mi_len - 1, max(set(length_count), key=length_count.count), ma_len + 1])
    ax.grid(axis='y')
    ax.set(title='Distribution of sequence lengths over all sequences', xlabel='Sequence Length(bp)')

    plt.savefig(imgname)

    # define report status
    status = 'good'
    if 0 in length_count:
        status = 'fail'
    elif min(length_count) != max(length_count):
        status = 'warning'

    return status


# -----------------------------------------------------------------------------
# create image and return status
def sequence_duplication_levels(sequences, fastq_name, imgname):
    bin_dict = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                7: 0, 8: 0, 9: 0, 10: 0, 50: 0, 100: 0,
                500: 0, 1000: 0, 5000: 0, 10000: 0}

    bin_dict_un = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0,
                   7: 0, 8: 0, 9: 0, 10: 0, 50: 0, 100: 0,
                   500: 0, 1000: 0, 5000: 0, 10000: 0}

    seq_count = {}
    
    def get_bin(val):
        if 10 <= val < 50:
            return 10
        elif 50 <= val < 100:
            return 50
        elif 100 <= val < 500:
            return 100
        elif 500 <= val < 1000:
            return 500
        elif 1000 <= val < 5000:
            return 1000
        elif 5000 <= val < 10000:
            return 5000
        elif 10000 <= val:
            return 10000
    
    def count_line(seq):
        if seq in seq_count:
            seq_count[seq] += 1
        else:
            seq_count[seq] = 1

    with open(fastq_name, 'r') as fastq_file:
        counter = 0
        reset_counter = False
        for line in fastq_file:
            counter += 1
            line = line[:50]
            if counter == 2 and not reset_counter:
                count_line(line)
                counter = 0
                reset_counter = True
            else:
                if counter % 4 == 0:
                    # line = line[:50]
                    count_line(line)
    for key, val in seq_count.items():
        if val < 10:
            bin_dict[val] += val
            bin_dict_un[val] += 1
        else:
            bin_dict[get_bin(val)] += val
            bin_dict_un[get_bin(val)] += 1

    x = ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10 ", ">50 ", ">100 ", ">500 ", ">1k ", ">5k ", ">10k "]
    y_total = [round((val / (counter // 4) + 1) * 100, 2) for val in bin_dict.values()]
    y_unique = [round((val / len(seq_count)) * 100, 2) for val in bin_dict_un.values()]
    # image creation example
    fig = plt.figure(figsize=(10, 7))
    # ax = fig.add_subplot(111)
    ax.set(title='sequence_duplication_levels')
    # removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    # make grid
    plt.grid(axis='x')
    plt.grid(axis='y')
    # ax.set_xlim([min_mq, max_mq])
    ax.set_ylim([0, 100])
    # add title
    ax.set(xlabel='Sequence Duplication Level', ylabel='')
    # add mean line
    ax.plot(x, y_total, color='blue')
    ax.plot(x, y_unique, color='red')
    fig.savefig(imgname)
    # define report status
    if y_total[1] > 50:
        status = 'fail'
    elif y_total[1] > 20:
        status = 'warning'
    else:
        status = 'good'
    return status


# -----------------------------------------------------------------------------
# create image and return status
def overrepresented_sequences(sequences, fastq_name, imgname):
    LINES_COUNT = 0
    SEQUENCES_SHA = {}
    SEQUENCES_PRINT = {}

    def process_line(line):
        seq_part = line[:50]
        line = line.encode('utf-8')
        line_sha = hashlib.sha256(line).hexdigest()
        if line_sha in SEQUENCES_SHA:
            if seq_part in SEQUENCES_PRINT:
                SEQUENCES_PRINT[seq_part] += 1
            else:
                SEQUENCES_PRINT[seq_part] = 2
        else:
            SEQUENCES_SHA[line_sha] = 1

    with open(fastq_name, 'r') as fastq_file:
        # read second line with sequence, after that process every 4th line
        reset_counter = False
        for line in fastq_file:
            LINES_COUNT += 1
            if LINES_COUNT == 2 and not (reset_counter):
                process_line(line)
                LINES_COUNT = 0
                reset_counter = True
            else:
                if LINES_COUNT % 4 == 0:
                    process_line(line)
    d_seqs = [k for k in SEQUENCES_PRINT.keys()]
    d_count = [v for v in SEQUENCES_PRINT.values()]
    data = {"Sequence": d_seqs,
            "Count": d_count}
    df = pd.DataFrame.from_dict(data)
    df["Percentage"] = round((df["Count"] / ((LINES_COUNT // 4) + 1)) * 100, 2)
    df = df.sort_values(by="Count", ascending=False)
    df_f = df[df["Percentage"] > 0.1]
    df_f = df_f.nlargest(n=20, columns=["Percentage"])
    fig = go.Figure(data=[go.Table(
        header=dict(values=list(df_f.columns),
                    fill_color='white',
                    align='left'),
        cells=dict(values=[df_f["Sequence"], df_f["Count"], df_f["Percentage"]],
                   fill_color='white',
                   align='left'))
    ])
    threshold = max(df_f["Percentage"])
    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    # define report status
    if 0.1 <= threshold < 1:
        status = 'warning'
    elif threshold > 1:
        status = 'fail'
    else:
        status = 'good'
    fig.update_layout(title="Top 20 of overrepresented sequences")
    fig.write_image(imgname)
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
