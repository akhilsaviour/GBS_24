import numpy as np
import quantp as qp
from functools import partial
import sax
##

def phot_cir(M=8,num_layer=8,wl=1.55):

    #num of inputs/output=m*2
    #width if mesh=n
    m=M//2
    n=num_layer//2
    t_no=m*n+(n-1)*(m)
    instances={}
    for t in range(1,t_no+1):
        var=f"{'S'}{t}"
        instances[var]='SU2'

    p_ar=[1,0]*m
    q_ar=[0,1]*m

    mesh_arr=[p_ar,q_ar]
    mesh=np.transpose(mesh_arr*n)
    mesh=mesh[:-1]   ## if the out put waveguide layer  should match with input layer
    # Mesh=mesh
    port_no=int(int(mesh.shape[0]+1)/2)
        
    conn_dict={}
    name_dict={}
    num=1
    for j in range(mesh.shape[1]):   # for j in range(mesh.shape[1]-1): 
        for i in range(mesh.shape[0]):
            if mesh[i][j]==1:
                string=f'({i},{j})'
                name_dict[str(num)]=[]
                name_dict[str(num)].append(string)
                num+=1
                if i==0:                    
                    conn_dict[string]=[]
                    conn_dict[string].append([[f'({i},{j+2})'],[f'({i+1},{j+1})']])
                elif i==mesh.shape[0]-1 :
                    
                    conn_dict[string]=[]
                    conn_dict[string].append([[f'({i-1},{j+1})'],[f'({i},{j+2})']])                    
                else:                    
                    conn_dict[string]=[]
                    conn_dict[string].append([[f'({i-1},{j+1})'],[f'({i+1},{j+1})']])

    conn={}
    port_num=int(int(mesh.shape[0]+1)/2)
    p1=1
    p2=port_num
    for ckey,cvalue in conn_dict.items():
        for nkey,nvalue in name_dict.items():      
            if ckey == nvalue[0]:
                # if int(nkey)>t_no-(port_num+(port_num-1)):             
                if int(nkey)==p1 :
                    string1=f"S{nkey},outt"
                    string2=f"S{nkey},outb"
                    conn[string1]=[]
                    conn[string2]=[]
                    p1+=(port_num+(port_num-1))                

                    for nkey1,nvalue1 in name_dict.items():
                        if cvalue[0][0]==nvalue1:
                            conn[string1]=f"S{nkey1},int"
                        if cvalue[0][1]==nvalue1:
                            conn[string2]=f"S{nkey1},int"

                elif  int(nkey)== p2:
                    string1=f"S{nkey},outt"
                    string2=f"S{nkey},outb"
                    conn[string1]=[]
                    conn[string2]=[]
                    
                    p2+=(port_num+(port_num-1))

                    for nkey1,nvalue1 in name_dict.items():
                        if cvalue[0][0]==nvalue1:
                            conn[string1]=f"S{nkey1},inb"
                        if cvalue[0][1]==nvalue1:
                            conn[string2]=f"S{nkey1},inb"     
                else:    
                    string1=f"S{nkey},outt"
                    string2=f"S{nkey},outb"
                    conn[string1]=[]
                    conn[string2]=[]

                    for nkey1,nvalue1 in name_dict.items():
                        if cvalue[0][0]==nvalue1:
                            conn[string1]=f"S{nkey1},inb"
                        if cvalue[0][1]==nvalue1:
                            conn[string2]=f"S{nkey1},int"
    pdict={}
    pind=1
    loop=2
    connections={}
    for kval,vval in conn.items():
        if vval==[]:
            if pind==loop:
                last_val=kval
                loop-=1
            else:
                pdict[pind]=[]
                pdict[pind]=[kval]
                pind+=1   
        else:
            connections[kval]=[]
            connections[kval]=[vval][0]

    pdict[pind]=[]
    pdict[pind]=[last_val]
    ports={}
    o=1
    for p in range(1,(port_no)+1):
        
        kin1=f"int{p}"
        kin2=f"inb{p}"

        kout1=f"outt{p}"
        kout2=f"outb{p}"

        ports[kin1]=[]
        ports[kin2]=[]
        ports[kin2]=[]
        ports[kout2]=[]

        ports[kin1]=f"S{p},int"
        ports[kin2]=f"S{p},inb"
        ports[kout1]=pdict[o][0]
        ports[kout2]=pdict[o+1][0]
        o+=2
    return instances,connections,ports



def sampler(a,b,c,wl,theta_arr=None,phi_arr=None):
       ##circuit
    if theta_arr is None  :
        theta_arr=np.zeros(len(a))
        print('theta array is None')
    if phi_arr is None:
        phi_arr=np.zeros(len(a))
        print('phi array is None')
    cir,info=sax.circuit(
        netlist={
            "instances":a,
            "connections":b,
            "ports":c
        },
        models={'SU2':partial(qp.trr)}) 
    
    su_dict={}
    i=0
    for key in range(len(a)):
        su_dict[key]={"theta":theta_arr[i],"phi":phi_arr[i]}
        i+=1
    
    val=cir(wl=wl,    S1=su_dict[0],
                      S2=su_dict[1],
                      S3=su_dict[2],
                      S4=su_dict[3],
                      S5=su_dict[4],
                      S6=su_dict[5],
                      S7=su_dict[6],
                      S8=su_dict[7],
                      S9=su_dict[8],
                      S10=su_dict[9],
                      S11=su_dict[10],
                      S12=su_dict[11],
                      S13=su_dict[12],
                      S14=su_dict[13],
                      S15=su_dict[14],
                      S16=su_dict[15],
                      S17=su_dict[16],
                      S18=su_dict[17],
                      S19=su_dict[18],
                      S20=su_dict[19],
                      S21=su_dict[20],
                      S22=su_dict[21],
                      S23=su_dict[22],
                      S24=su_dict[23],
                      S25=su_dict[24],
                      S26=su_dict[25],
                      S27=su_dict[26],
                      S28=su_dict[27],
                      
    )
    return val
