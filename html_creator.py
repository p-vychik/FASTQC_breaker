
class Html:
    ''' Creator of HTML page '''

    def __init__(self, template_filename):
        # read template file
        with open(template_filename, 'r') as f:
            self.filedata = f.read()

    def replace_words(self, words):
        # replace the target string
        for key in words:
            # print(f'  Replace:  {key} -> {words[key]}')
            self.filedata = self.filedata.replace(key, words[key])

    def make_report(self, report_name):
        # write the report file
        with open(report_name, 'w') as f:
            f.write(self.filedata)
