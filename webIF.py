"""
WebIF is a class that handles sending commands to the 4D Web Service that enables
remote measurements of C1S04. In order for WebIF to work properly, ensure:
1.) The interferometer is turned on.
2.) 4Sight is active.
3.) The optic is aligned, and has a mask that's been set manually.
4.) In 4Sight click: Tools->Start Web Service
"""

import requests

class WebIF:

    def __init__(self):
        self.url = 'http://localhost/WebService4D/WebService4D.asmx'
        self.timeout = 10000 # default timeout value of Web Service

    def setTimeout(self, t):
        """
        Set the amount of time to wait for response from server, in ms.
        t: integer
        """
        response = requests.get(self.url+'/SetTimeout?timeOut={}'.format(int(t)))
        returnStr = response.text[-16:-9]
        if returnStr == 'success': self.timeout = t
        else: returnStr = 'fail'
        print('setTimeout result: {}'.format(returnStr))

    def getTimeout(self):
        """
        Returns the current timeout value, in ms, used by the Web Service.
        """
        response = requests.get(self.url+'/GetTimeout?')
        trimmed_response = response.text.split('>')[-2]
        end_idx = trimmed_response.index('<')
        timeVal = int(trimmed_response[:end_idx])
        return timeVal

    def averageMeasure(self, count, printout=False):
        """
        Completes an average measurement using a number of measurements
        specified by count.
        count: integer
        """
        response = requests.get(self.url+'/AverageMeasure?count={}'.format(int(count)))
        # print(response.text)
        returnStr = response.text[-16:-9]
        if returnStr != 'success': returnStr = 'fail'
        if printout:
            print('averageMeasure result: {}'.format(returnStr))

    def saveMeasurement(self, fname, printout=False):
        """
        Saves the current measurement as a .h5 file. fname should be the full
        path for the file to be saved using '\\' Python-style delimiters between
        directories.
        """
        response = requests.get(self.url+'/SaveMeasurement?fileName={}'\
                                .format(fname.replace("\\", "/")))
        returnStr = response.text[-16:-9]
        if returnStr != 'success': returnStr = 'fail'
        if printout:
            print('saveMeasurement result: {}'.format(returnStr))
