import time
from collections import defaultdict

class timer:
    def __init__(self, running_log, rank, size):
        self.running_log = running_log
        self.rank = rank
        self.size = size
        self.task_time = defaultdict(dict)
        self.program_start_time = None

    def program_start(self):
        ttime = time.asctime(time.localtime(time.time()))
        if self.rank == 0:
            with open(self.running_log, 'w') as f:
                print('program start at', ttime, file=f)
                print('processor number is', self.size, file=f)

        self.program_start_time = time.time()

    def program_end(self):
        total_time = time.time() - self.program_start_time
        if self.rank == 0:
            with open(self.running_log, 'a') as f:
                h, m = divmod(total_time, 3600)
                m, s = divmod(m, 60)
                s = int(s)
                print('Total time', ': %d h %d m %d s' % (h, m, s), file=f)

    def moment(self, description):
        ttime = time.asctime(time.localtime(time.time()))
        if self.rank == 0:
            with open(self.running_log, 'w') as f:
                print(description, ttime, file=f)

    def start(self, class_name, description):
        start_time = time.time()
        temp = {description : [start_time]}
        self.task_time[class_name].update(temp)


    def end(self, class_name, description):
        start_time = self.task_time[class_name][description][0]
        end_time = time.time()
        duration = end_time - start_time
        self.task_time[class_name][description].extend([end_time, duration])

    def print_all(self):
        if self.rank == 0:
            maxlen = 0
            for description in self.task_time.values():
                tem_maxlen = max(map(len, description.keys()))
                maxlen = max(maxlen, tem_maxlen)
            with open(self.running_log, 'a') as f:
                f.write('\n')
                f.write('\n------------------------------------------------------')
                f.write('\n|                                                    |')
                f.write('\n|                   Time Statistics                  |')
                f.write('\n|                                                    |')
                f.write('\n------------------------------------------------------')
                f.write('\n\n')
                for class_name, description in self.task_time.items():
                    print(class_name + ' : ', file=f)
                    for key, value in description.items():
                        print(' >> ' + key.ljust(maxlen), ': %15.6e s' % value[2], file=f)
