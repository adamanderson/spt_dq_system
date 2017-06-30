import sys
from os import listdir

def main():
    # load arguments into a dict to be passed to plotting functions
    request = {'type': sys.argv[1], 'path': sys.argv[2], 'source': sys.argv[3], 'observation': sys.argv[4]}

    # load python file and function
    module = __import__(request['type'])
    func = getattr(module, request['type'] + '_plot')
    # call plotting function print plot name without any newline
    func(request)
    sys.stdout.write('plot.png')

#start process
if __name__ == '__main__':
    main()
