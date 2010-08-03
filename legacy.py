import types

import vdj

# def parse_VDJXML(inputfile):
#     """Generator to return ImmuneChain objects from a vdjxml file.
#     
#     Utilizes python XML libraries.
#     """
#     # NOTE: this can probably be made more elegant if implemented
#     # as a class
#     
#     # global state variables.
#     chain = None
#     data_buffer = None
#     
#     def start_handler(name,attributes):
#         if name == 'ImmuneChain':
#             chain = ImmuneChain()
#         else:
#             data_buffer = ''
#     
#     def end_handler(name):
#         if name == 'ImmuneChain':
#             yield chain
#         elif name == 'tag':
#             chain.add_tags(data_buffer)
#         else:
#             chain.__setattr__(name,data_buffer)
#     
#     def data_handler(data):
#         data_buffer += data
#     
#     xmlparser = xml.parsers.expat.ParserCreate()
#     xmlparser.StartElementHandler  = start_handler
#     xmlparser.EndElementHandler    = end_handler
#     xmlparser.CharacterDataHandler = data_handler
#     
#     if not hasattr(inputfile,'read'):
#         inputfile = open(inputfile,'r')
#     
#     xmlparser.ParseFile(inputfile)


# =============================================================
# = START DEPRECATED ==========================================
# =============================================================

def parse_VDJXML_old(inputfile):
    """Load a data from a VDJXML file as a Repertoire or list of ImmuneChains
    
    NOTE: this fn does NOT utilize the XML libraries; it implements a manual parser
    that takes input line by line.
    
    THIS ASSUMES THAT EVERY XML ELEMENT TAKES ONE AND ONLY ONE LINE
    
    """
    
    if isinstance(inputfile,types.StringTypes):
        ip = open(inputfile,'r')
    elif isinstance(inputfile,file):
        ip = inputfile
    
    numChains = 0
    
    possible_elements = [
                'descr',
                'seq',
                'v',
                'd',
                'j',
                'ighc',
                'cdr3',
                'junction',
                'func',
                'tag'
                ]
    
    for line in ip:
        line = line.strip()
        endelementpos = line.find('>') + 1
        xmlelement = line[0:endelementpos]
        element = xmlelement[1:-1]
        
        if xmlelement == '<ImmuneChain>':
            chain = vdj.ImmuneChain()
        elif xmlelement == '</ImmuneChain>':
            numChains += 1
            yield chain
        elif element in possible_elements:
            if element == 'tag':
                tagdata = line[endelementpos:-1*(endelementpos+1)].split('|')
                if tagdata[0] == 'experiment':
                    chain.__setattr__('experiment','|'.join(tagdata[1:]))
                elif tagdata[0] == 'clone':
                    chain.__setattr__('clone','|'.join(tagdata[1:]))
                elif tagdata[0] == 'barcode':
                    chain.__setattr__('barcode','|'.join(tagdata[1:]))
                elif tagdata[0] == 'v_end_idx':
                    chain.__setattr__('v_end_idx','|'.join(tagdata[1:]))
                elif tagdata[0] == 'j_start_idx':
                    chain.__setattr__('j_start_idx','|'.join(tagdata[1:]))
                else:
                    chain.add_tags(line[endelementpos:-1*(endelementpos+1)])
            elif element == 'ighc':
                eltdata = line[endelementpos:-1*(endelementpos+1)]
                if eltdata != '':
                    chain.c = eltdata
                else:
                    pass
            elif element == 'func':
                pass
            else:
                eltdata = line[endelementpos:-1*(endelementpos+1)]
                if eltdata != '':
                    chain.__setattr__(element,eltdata)
                else:
                    pass
    
    if isinstance(inputfile,types.StringTypes):
        ip.close()

# =============================================================
# = END DEPRECATED ============================================
# =============================================================



# def filter_parse_VDJXML(inputfile,predicate):
#     """Load a data from a VDJXML file as a Repertoire or list of ImmuneChains
#     
#     predicate is a function that takes a chain and return True or False.  Things
#     that return false are skipped.
#     
#     NOTE: this fn does NOT utilize the XML libraries; it implements a manual parser
#     that takes input line by line.
#     
#     THIS ASSUMES THAT EVERY XML ELEMENT TAKES ONE AND ONLY ONE LINE
#     
#     """
#     
#     if isinstance(inputfile,types.StringTypes):
#         ip = open(inputfile,'r')
#     elif isinstance(inputfile,file):
#         ip = inputfile
#     
#     numChains = 0
#     
#     possible_elements = [
#                 'descr',
#                 'seq',
#                 'v',
#                 'd',
#                 'j',
#                 'ighc',
#                 'cdr3',
#                 'junction',
#                 'func',
#                 'tag'
#                 ]
#     
#     for line in ip:
#         line = line.strip()
#         endelementpos = line.find('>') + 1
#         xmlelement = line[0:endelementpos]
#         element = xmlelement[1:-1]
#         
#         if xmlelement == '<ImmuneChain>':
#             chain = ImmuneChain()
#         elif xmlelement == '</ImmuneChain>':
#             numChains += 1
#             if predicate(chain) == True:
#                 yield chain
#             else:
#                 pass
#         elif element in possible_elements:
#             if element == 'cdr3':
#                 chain.cdr3 = eval(line[endelementpos:-1*(endelementpos+1)])
#             elif element == 'tag':
#                 chain.add_tags(line[endelementpos:-1*(endelementpos+1)])
#             else:
#                 chain.__setattr__(element,line[endelementpos:-1*(endelementpos+1)])
#     
#     if isinstance(inputfile,types.StringTypes):
#         ip.close()

