import shutil
import time
import sys
import os
import sequences as sq
import html_creator
import reports


# -----------------------------------------------------------------------------
def do_check(input_fastq, output_dir):

    tm1 = time.time()

    # open files
    fi = open(input_fastq, "r")

    # parse line by line
    lines = fi.readlines()

    # find max sequence length
    read_num = int(len(lines) / 4)
    max_seq_len = 0
    for i in range(read_num):
        l2 = lines[i * 4 + 1].strip()
        max_seq_len = max(max_seq_len, len(l2))
    # print(f'max_seq_len:  {max_seq_len}')

    # allocate Sequences for 'read_num' and 'max_seq_len'
    sequences = sq.Sequences(max_seq_len, read_num)

    print('Parse file:')
    point_dot = 0
    point_percent = -1
    for i in range(read_num):
        l1 = lines[i * 4 + 0].strip()
        l2 = lines[i * 4 + 1].strip()
        l4 = lines[i * 4 + 3].strip()
        sequences.add(l1, l2, l4)
        # print progress
        new_dot = int(i * 50 / read_num)
        new_percent = int(i * 10 / read_num)
        if new_percent != point_percent:
            point_percent = new_percent
            sys.stdout.write(str(point_percent*10) + '%')
            sys.stdout.flush()
        elif new_dot != point_dot:
            point_dot = new_dot
            sys.stdout.write('.')
            sys.stdout.flush()
    print('100%')

    # sequences.print()
    tm2 = time.time()
    print(f'1.  parse file:                 {(tm2 - tm1):.3f} sec')

    # Make reports

    html = html_creator.Html('./Html/template.html')
    html.replace_words({'FASTQ_FILE_NAME': input_fastq})

    # Basic Statistics
    words = reports.basic_statistics(sequences, input_fastq)
    tm3 = time.time()
    print(f'2.  basic_statistics:           {(tm3 - tm2):.3f} sec')
    html.replace_words(words)

    # Per base sequence quality
    img = output_dir + '/per_base_sequence_quality.png'
    status = reports.per_base_sequence_quality(sequences, input_fastq, img)
    tm4 = time.time()
    print(f'3.  per_base_seq_quality:       {(tm4 - tm3):.3f} sec')
    html.replace_words({'PER_BASE_SEQUENCE_QUALITY_status': status})

    # Per sequence quality scores
    img = output_dir + '/per_sequence_quality_scores.png'
    status = reports.per_sequence_quality_scores(sequences, input_fastq, img)
    tm5 = time.time()
    print(f'4.  per_seq_quality_scores:     {(tm5 - tm4):.3f} sec')
    html.replace_words({'PER_SEQUENCE_QUALITY_SCORES_status': status})

    # Per base sequence content
    img = output_dir + '/per_base_sequence_content.png'
    status = reports.per_base_sequence_content(sequences, input_fastq, img)
    tm6 = time.time()
    print(f'5.  per_base_seq_content:       {(tm6 - tm5):.3f} sec')
    html.replace_words({'PER_BASE_SEQUENCE_CONTENT_status': status})

    # Per base GC content
    img = output_dir + '/per_base_gc_content.png'
    status = reports.per_base_gc_content(sequences, input_fastq, img)
    tm6 = time.time()
    print(f'6.  per_base_gc_content:        {(tm6 - tm5):.3f} sec')
    html.replace_words({'PER_BASE_GC_CONTENT_status': status})

    # Per sequence GC content
    img = output_dir + '/per_sequence_gc_content.png'
    status = reports.per_sequence_gc_content(sequences, input_fastq, img)
    tm7 = time.time()
    print(f'7.  per_sequence_gc_content:    {(tm7 - tm6):.3f} sec')
    html.replace_words({'PER_SEQUENCE_GC_CONTENT_status': status})

    # Per base N content
    img = output_dir + '/per_base_n_content.png'
    status = reports.per_base_n_content(sequences, input_fastq, img)
    tm8 = time.time()
    print(f'8.  per_base_n_content:         {(tm8 - tm7):.3f} sec')
    html.replace_words({'PER_BASE_N_CONTENT_status': status})

    # Sequence Length Distribution
    img = output_dir + '/sequence_length_distribution.png'
    status = reports.sequence_length_distribution(sequences, input_fastq, img)
    tm9 = time.time()
    print(f'9.  seq_length_distribution:    {(tm9 - tm8):.3f} sec')
    html.replace_words({'SEQUENCE_LENGTH_DISTRIBUTION_status': status})

    # Sequence Duplication Levels
    img = output_dir + '/sequence_duplication_levels.png'
    status = reports.sequence_duplication_levels(sequences, input_fastq, img)
    tm10 = time.time()
    print(f'10.  seq_duplication_levels:    {(tm10 - tm9):.3f} sec')
    html.replace_words({'SEQUENCE_DUPLICATION_LEVELS_status': status})

    # Overrepresented sequences
    img = output_dir + '/overrepresented_sequences.png'
    status = reports.overrepresented_sequences(sequences, input_fastq, img)
    tm11 = time.time()
    print(f'11. overrepresented_sequences:  {(tm11 - tm10):.3f} sec')
    html.replace_words({'OVERREPRESENTED_SEQUENCES_status': status})

    # Adapter Content
    img = output_dir + '/adapter_content.png'
    status = reports.adapter_content(sequences, input_fastq, img)
    tm12 = time.time()
    print(f'12. adapter_content:            {(tm12 - tm11):.3f} sec')
    html.replace_words({'ADAPTER_CONTENT_status': status})

    # generate html report
    rep_name = output_dir + '/report.html'
    html.make_report(rep_name)
    print(f'13. report created:             {rep_name}')

    # close files
    fi.close()


def get_script_parameters():

    input_fastq = ''
    output_dir = ''

    if len(sys.argv) < 5:
        print(f'USE:  {sys.argv[0]} -i <fastq-file-name> -o <dir-for-report>')
        sys.exit(-1)  # exit programm with error code

    if sys.argv[1] == '-i' or sys.argv[1] == '--input':
        input_fastq = sys.argv[2]

    if sys.argv[3] == '-i' or sys.argv[3] == '--input':
        input_fastq = sys.argv[4]

    if sys.argv[1] == '-o' or sys.argv[1] == '--output':
        output_dir = sys.argv[2]

    if sys.argv[3] == '-o' or sys.argv[3] == '--output':
        output_dir = sys.argv[4]

    if input_fastq == '' or output_dir == '':
        print(f'USE:  {sys.argv[0]} -i <fastq-file-name> -o <dir-for-report>')
        sys.exit(-1)  # exit programm with error code

    return input_fastq, output_dir


# -----------------------------------------------------------------------------
# main function with default values
def main():

    input_fastq, output_dir = get_script_parameters()

    try:
        shutil.rmtree(output_dir)  # del report dir with content
    except FileNotFoundError:
        print(end='')              # ignore, do nothing
    os.mkdir(output_dir)           # create new report dir
    shutil.copyfile('Html/good.png',    output_dir + '/good.png')
    shutil.copyfile('Html/warning.png', output_dir + '/warning.png')
    shutil.copyfile('Html/fail.png',    output_dir + '/fail.png')

    # read file and create reports as html
    do_check(input_fastq, output_dir)


# execute main function
main()
