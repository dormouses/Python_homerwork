#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import Tkinter as tk
import Verle as v

import time
import random as rand
import sys

def add():
    m=scale_m.get()
    x=text_x.get()
    y=text_y.get()
    vx=text_vx.get()
    vy=text_vy.get()
    t=text_t.get()
    c=rand.randint(0, len(collors) - 1)
    system.make_part(int(x), int(y), int(vx), int(vy), int(m), int(t), collors[c])
    #print 'hey'

def show():
        s=0
        c=0
        while True:
            time.sleep(0.5)
            picture.delete("all")
            for i in system.parts:
                picture.create_oval(i.coord.x-i.r, i.coord.y-i.r, i.coord.x+i.r, i.coord.y+i.r, fill=i.col)   
            picture.update()
            start = time.time()
            system.update() 
            end = time.time()
            c+=1
            s+=end-start
            if c==100:
                print s/c
                c=0
            #print end-start
            #print system.type
            sys_type=system_type.get()
            if sys_type=="Verle":
                system.type=0
            elif sys_type=="Parallel":
                system.type=1 
            elif sys_type=="Odeint":
                system.type=2
            elif sys_type=="Cython":
                system.type=3
            
            
        
sys.path.insert(0, '/Users/maria/Desktop/pyt/5_ball')
  
collors=['#800000', '#ff0000', '#ff00ff', '#00ffff', '#000000', '#00ff00', '#000080', '#ffff00'] 

root = tk.Tk()
root.title("Шарики")
root.geometry('800x450')
  
space = 10
w_obj = 35
line = 450
w_lable = 100
w_text = 200
picture = tk.Canvas(width=400, height = 440)
picture.place(x=5, y=5)

c_y=space
label_x = tk.Label(text="X")
label_x.place(x=line, y=c_y, width=w_lable, height=w_obj)
text_x = tk.Entry()
text_x.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
label_y = tk.Label(text="Y")
label_y.place(x=line, y=c_y, width=w_lable, height=w_obj)
text_y = tk.Entry()
text_y.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
label_vx = tk.Label(text="Vx")
label_vx.place(x=line, y=c_y, width=w_lable, height=w_obj)
text_vx = tk.Entry()
text_vx.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
label_vy = tk.Label(text="Vy")
label_vy.place(x=line, y=c_y, width=w_lable, height=w_obj)
text_vy = tk.Entry()
text_vy.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
label_t = tk.Label(text="Time")
label_t.place(x=line, y=c_y, width=w_lable, height=w_obj)
text_t = tk.Entry()
text_t.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
label_m = tk.Label(text="Wigth")
label_m.place(x=line, y=c_y, width=w_lable, height=w_obj)
scale_m = tk.Scale(orient="horizontal", from_=1000,to=2000,tickinterval=50, resolution=50)
scale_m.place(x=line+w_lable, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
all_menu = ["Verle", "Odeint", "Parallel", "Cython"]
system_type = tk.StringVar()
system_type.set("Verle")
menu = tk.OptionMenu(root, system_type, *all_menu)
menu.place(x=line+w_lable//2, y=c_y, width=w_text, height=w_obj)

c_y+=w_obj+space+space
but_add = tk.Button(text="Добавить", command=add)
but_add.place(x=line+w_lable//2, y=c_y, width=w_text, height=w_obj)

system=v.Verle()
root.after(0, show)
root.mainloop()
