import numpy as np
import matplotlib.pyplot as plt
import hashlib
import plotly.graph_objects as go
import pandas as pd


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
def sequence_length_distribution(sequences, fastq_name, imgname):

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
    #***
    # image creation example
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    ax.set(title='sequence_duplication_levels')
    # removing top axes and right axes ticks
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # make grid
    plt.grid(axis='x')
    plt.grid(axis='y')
    #ax.set_xlim([min_mq, max_mq])
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
