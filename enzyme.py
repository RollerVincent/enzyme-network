import requests
import os


class Enzymes:
    data = None
    counts = {}
    reactands = {}
    exclusion = []
    def __init__(self):
        self.download()
        self.parse()
        self.filter()
        self.prepare()
        self.stats()
      
    def download(self):
        if not os.path.exists('./data'):
            os.mkdir('./data')
        if not os.path.exists('./data/enzyme.dat'):
            request = requests.get('https://ftp.expasy.org/databases/enzyme/enzyme.dat', stream=True)
            with open('./data/enzyme.dat', 'w') as w:
                for l in request.iter_lines():
                    w.write(l.decode() + '\n')

    def parse(self):
        with open('./data/enzyme.dat') as f:
            ec = None
            ecs = []
            d = {}
            for l in f:
                if l[0] == '/':
                    if ec is not None:
                        ec.data = d
                        ecs.append(ec)
                    ec = EC()
                    d = {}
                else:
                    if ec is not None:
                        c = l[:2]
                        if c not in d:
                            d.update({c:''})
                        d[c] += l[5:-1]
            ec.data = d
            ecs.append(ec)
        self.data = ecs
        print(str(len(self.data)) + ' enzyme classes')

    def filter(self):
        tmp = []
        for ec in self.data:
            if 'CA' in ec.data:
                ca = ec.data['CA'][:-1].split(' = ') #   handle ('<=>')
                if(len(ca) == 2):
                    ec.left = ca[0]
                    ec.right = ca[1]
                    tmp.append(ec)
        self.data = tmp
        print(str(len(self.data)) + ' enzyme classes with clear reactions')

    def prepare(self):
        for ec in self.data:
            for i in range(2):
                d = None
                if i==0:
                    d = ec.left
                else:
                    d = ec.right
                
                status = True
                for j,letter in enumerate(d):
                    if letter == '(':
                        status = False
                    elif letter == ')':
                        status = True
                    elif letter == '+' and status == True:
                        d = d[:j] + '$' + d[j+1:]
                
                d = [x.strip().upper() for x in d.split('$')]
                for r in d:
                    n = r.split(' ')[0]
                    if n.isdigit():
                        r = r[len(n)+1:]

                tmp = set()
                for r in d:
                    n = r.split(' ')[0]
                    if n.isdigit() or n in ['AN', 'A']:
                        tmp.add(r[len(n)+1:])
                    else:
                        tmp.add(r)
                d = tmp
                if i == 0:
                    ec.left = d
                else:
                    ec.right = d

    def stats(self):
        for ec in self.data:
            for i,d in enumerate([ec.left, ec.right]):
                for r in d:
                    if r not in self.counts:
                        self.counts.update({r:[0,0]})
                        self.reactands.update({r:[set(),set()]})
                    self.counts[r][i] += 1
                    self.reactands[r][i].add(ec)
        self.counts = sorted(self.counts.items(), key=lambda kv: -(kv[1][0] + kv[1][1]))

    def exclude(self):
        tmp = {}
        for r in self.reactands:
            if r not in self.exclusion:
                tmp.update({r:self.reactands[r]})
        self.reactands = tmp

    def build(self):
        self.exclude()
        c = 0
        for r in self.reactands:
            reaction = self.reactands[r]
            for enters in reaction[0]:
                for exits in reaction[1]:
                    exits.children.update({r:enters})
                    enters.parents.update({r:exits})
                    c += 1

        print('graph with '+str(len(self.data))+' nodes and '+str(c)+' links.')

    def save(self):
        with open('enzyme_network.gnf', 'w') as w:
            w.write('NODE\tID\tLABEL\n')
            for ec in self.data:
                w.write('NODE\t' + ec.data['ID'] + '\t' + ec.data['DE'][:-1] + '\n')
            w.write('\nLINK\tSRC\tDST\tLABEL\n')
            for ec in self.data:
                for r in ec.children:
                    ch = ec.children[r]
                    w.write('LINK\t' + ec.data['ID'] + '\t' + ch.data['ID'] + '\t' + r + '\n')

class EC:
    data = None
    left = None
    right = None
    parents = {}
    children = {}


enzymes = Enzymes()
enzymes.exclusion = [
                     'H(2)O', 
                     'NADP(+)', 
                     'H(+)', 
                     'NADPH', 
                     'ATP', 
                     'ADP', 
                     'DIPHOSPHATE', 
                     'PHOSPHATE', 
                     'NAD(+)', 
                     'NADH', 
                     'CO(2)'
                    ]
enzymes.build()
enzymes.save()