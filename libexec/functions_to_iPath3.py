import sys
import requests
from IPython.display import SVG
import svgutils.compose as sc


def add_color(value):
    green_gradient = ["#d8fcd8", "#a0f9a0", "#62f163", "#2aed2d", "#03df04"]
    red_gradient = ["#F8D1CD", "#F0A8AB", "#E97E88", "#E15566", "#DA2C43"]
    if value < 0:
        gradient = green_gradient
    else:
        gradient = red_gradient
    if value < 3:
        return gradient[0]
    elif value >= 3 and value < 6:
        return gradient[1]
    elif value >= 6 and value < 9:
        return gradient[2]   
    elif value >= 9 and value < 12:
        return gradient[3]   
    elif value >= 12:
        return gradient[4]

def add_width(value):
    if value < 0.00001:
        return "W10"
    elif value < 0.0001:
        return "W8"
    elif value < 0.001:
        return "W6"
    elif value < 0.01:
        return "W4"
    else:
        return 'W2'

def to_parameters(selection ,
                          export_type='svg',
                          include_metabolic=True,
                          include_secondary=False,
                          include_antibiotic=False,
                          include_microbial=False,
                          whole_modules=False,
                          whole_pathways=False,
                          keep_colors=False,
                          default_opacity=1,
                          default_width=3,
                          default_radius=7,
                          default_color='#666666',
                          query_reactions=False,
                          tax_filter='',
                          export_dpi=1200):

    allowed_export_types= ['svg','png','pdf','eps']
    assert export_type in allowed_export_types , "export_type {} needs to be one of {}".froamt(export_type,allowed_export_types)
    #assert map_type=='svg', "I can not save PNG images"

    return dict( selection=selection,
                export_type=export_type,
                keep_colors=int(keep_colors),
                include_metabolic=int(include_metabolic),
                include_secondary= int(include_secondary),
               include_antibiotic= int(include_antibiotic),
               include_microbial= int(include_microbial),
               whole_modules= int(whole_modules),
               default_opacity= default_opacity,
                whole_pathways=int(whole_pathways),
               default_width= default_width,
               default_color= default_color,
               default_radius= default_radius,
               query_reactions= int(query_reactions),
               tax_filter= tax_filter,
               export_dpi=export_dpi)

def get_map(selection,map_name='map',**param):

    url= 'https://pathways.embl.de/mapping.cgi'



    #print(selection[:300]+'...')
    parameters=to_parameters(selection,**param)

    r = requests.post(url, data=parameters)
    assert r.ok, r.text


    with open(map_name+'.svg','w') as file:
        file.write(r.text)


def scale_map(map_name):

    sc.Figure("26cm", "16cm",
        sc.Panel(sc.SVG(map_name+".svg").scale(0.265)),
        ).save(map_name+'_scaled.svg')
        
def generate_ipath_input( FH_in ):
    out_str = str()
    with open("../../tools/deseq2_visualisation/over.tsv") as FH_in:
        for li in FH_in:
            if li.startswith("OTU"):
                continue
            else:
                li = li.strip().split('\t')
                func = li[0]
                l2fc = float(li[1])
                padj = float(li[2])
                out_str += ' '.join([ li[0], add_color(l2fc), add_width(padj) ]) + '\n'
    return out_str

str = generate_ipath_input( "../../tools/deseq2_visualisation/over.tsv")
get_map(str)

scale_map('map')
