from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
from tkinter import *
from PIL import Image
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import tornado

root = tk.Tk()
root.geometry("1050x680")
root.resizable(0,0)

# ======================================================================#
# 						EXTERNAL FUNCTION 								#
# ======================================================================#

def initial_Condition(x, A0, xC, xwide) :
	return A0*np.exp((-((x-xC)/xwide)*((x-xC)/xwide)))

def call_SWE(xleft, xright, Nx, dx, tfin, A0, xC, xwide, hC, h0, hshall, hwide) :
	# 0. INPUT PARAMETERS
	# wave initial condition
	eta_list = list()
	time_list = list()
	x  = [0]*(Nx+1)
	x  = np.linspace(xleft, xright, Nx + 1)

	dt      = 0.0001
	# parameters
	g       = 9.81

	# 1. PREPARATION 
	h  		= [0]*(Nx+1)
	H  		= [0]*(Nx+1)	
	H_star  = [0]*(Nx+2)  # total depth in 'half-grid'
	H_bar  	= [0]*(Nx+2)  # total depth at half grid
	u 		= [0]*(Nx+2)  # velocity in 'half-grid'
	u_star  = [0]*(Nx+1)
	u_du 	= [0]*(Nx+2)  # advection term in half-grid
	u_new  	= [0]*(Nx+2)  # updated velocity u
	q  		= [0]*(Nx+2)  # q flux = H*u at half grid
	q_bar  	= [0]*(Nx+1)  # q at full grid
	eta  	= [0]*(Nx+1)  # old eta
	eta_new = [0]*(Nx+1)  # updated eta

	eta = initial_Condition(x, A0, xC, xwide)
	h = h0 - hshall*np.exp(-((x-hC)/hwide)*((x-hC)/hwide)) # 2.2. Bathymetry/bottom
	H = h + eta                                            # 2.3. total depth in full grid

	t = 0
	step = 0

	while (t <= tfin) :
		t = t + dt
		step += 1
		for j in range(1, Nx+1) :
			if (u[j] >= 0) :
				H_star[j]     = H[j-1]
			else :
				H_star[j]     = H[j]
			q[j] = H_star[j]*u[j]
		H_star[0]       = H_star[1]
		H_star[Nx+1]    = H_star[Nx]
		q[0]    = 0
		q[Nx+1] = 0


		for j in range(0, Nx+1) :
			if(j==0 or j==Nx) :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
			else :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
			H[j] = eta_new[j] + h[j]
			q_bar[j] = 0.5*(q[j+1] + q[j])
			if(q_bar[j] >= 0) :
				u_star[j]=u[j]
			else :
				u_star[j]=u[j+1]
		# 3.2. Calculate 2nd eq. (momentum eq.)
		for j in range(1, Nx+1) :
			H_bar[j] = 0.5*(H[j] + H[j-1])
			if(H_bar[j] == 0) :
				u_new[j] = 0
			else :
				u_du[j] = 1./(H_bar[j]*dx)*((q_bar[j]*u_star[j] - q_bar[j-1]*u_star[j-1]) - (u[j]*(q_bar[j] - q_bar[j-1])))
				u_new[j] = u[j] - (dt/dx*g*( eta_new[j] - eta_new[j-1] )) - (dt*u_du[j])
		u_new[0] = 0
		u_new[Nx+1] = 0

		# updating values of H, eta, u
		eta   = eta_new
		H     = h + eta
		u     = u_new
		if (step % 2 == 0) :
			eta_list.append(eta_new[:])
			time_list.append(t)
	return eta_list

# ======================================================================#
# 						EXTERNAL FUNCTION 								#
# ======================================================================#

def on_closing():
    if messagebox.askokcancel("Quit", "Do you want to Quit?"):
        root.destroy()

class Window(tk.Frame):

	def __init__(self, master=None):
		tk.Frame.__init__(self, master)
		self.master = master
		self.init_window()
		self.running = False
		self.animate = None

	def start(self):
		self.eta_list = list()
		self.time_list = list()

		self.Axleft=float(self.xleft.get())
		self.Axright=float(self.xright.get())
		self.ANx=int(self.Nx.get())
		self.Adx=float(self.dx.get())
		self.Atfin=float(self.tfin.get())
		self.AA0=float(self.A0.get())
		self.AxC=float(self.xC.get())
		self.Axwide=float(self.xwide.get())
		self.Adt=float(self.dt.get())
		self.AhC=float(self.hC.get())
		self.Ah0=float(self.h0.get())
		self.Ahshall=float(self.hshall.get())
		self.Ahwide=float(self.hwide.get())

		self.x  = [0]*(self.ANx+1)
		self.x = np.linspace(self.Axleft, self.Axright, self.ANx + 1)
		self.eta_list = call_SWE(self.Axleft, self.Axright, self.ANx, self.Adx, self.Atfin, self.AA0, self.AxC, self.Axwide, self.AhC, self.Ah0, self.Ahshall, self.Ahwide)
		self.eta_len = len(self.eta_list)

		def animate(i):
			self.y = self.eta_list[i][:]
			self.line.set_data(self.x,self.y)  # update the data
			p = plt.fill_between(self.x, -1, self.y, facecolor = 'C0', alpha = 0.2)
			return self.line, p

		self.fig = plt.figure()
		self.ax = self.fig.add_subplot(1,1,1)
		self.ax.set_xlim([0,1])
		self.ax.set_ylim([-0.002, 0.002])
		self.line, = self.ax.plot(np.add([], []), lw=5, color='black')
		self.ax.legend(['Analytics'])
		self.ax = plt.grid()
		plt.xlabel("x (m)",fontsize=14)
		self.ax = plt.title("\nMOmentum conservation Saint-Venant 1 Dimension")
		self.fig.suptitle('Simulation MoSV 1D', fontsize=20)

		self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame6)
		self.canvas.draw()
		self.canvas.get_tk_widget().grid(column=3,row=1, padx=10)

		self.ani = animation.FuncAnimation(self.fig, animate, frames=np.arange(1, self.eta_len, 10), interval=30, blit=True, cache_frame_data=False, repeat_delay=1000)
		self.running = True
		self.B5.config(text='Pause')
		self.ani._start()
	
	def on_click(self):

		if self.ani is None:
			return self.start()
		if self.running:
			self.ani.event_source.stop()
			self.B5.config(text='Run')
		else:
			self.ani.event_source.start()
			self.B5.config(text='Pause')
		self.running = not self.running

	def clear_text(self):
		self.E0.delete(0, 'end')
		self.E1.delete(0, 'end')
		self.E2.delete(0, 'end')
		self.E3.delete(0, 'end')
		self.E4.delete(0, 'end')
		self.E5.delete(0, 'end')
		self.E6.delete(0, 'end')
		self.E7.delete(0, 'end')
		self.E8.delete(0, 'end')
		self.E9.delete(0, 'end')
		self.E10.delete(0, 'end')
		self.E11.delete(0, 'end')
		self.E12.delete(0, 'end')

	def validate_float(self, action, index, value_if_allowed,
	    prior_value, text, validation_type, trigger_type, widget_name):
	    # action=1 -> insert
	    if(action=='1'):
	        if text in value_if_allowed:
	            try: 
	                float(value_if_allowed)
	                return True
	            except ValueError:
	                return False
	        else:
	            return False
	    else:
	        return True

	def save(self):
	    file = self.Ename.get()
	    param0=self.E0.get()
	    param1=self.E1.get()
	    param2=self.E2.get()
	    param3=self.E3.get()
	    param4=self.E4.get()
	    param5=self.E5.get()
	    param6=self.E6.get()
	    param7=self.E7.get()
	    param8=self.E8.get()
	    param9=self.E9.get()
	    param10=self.E10.get()
	    param11=self.E11.get()
	    param12=self.E12.get()
	    with open(file + '.txt', 'w') as file_object:
	        file_object.write(param0)
	        file_object.write(',')
	        file_object.write(param1)
	        file_object.write(',')
	        file_object.write(param2)
	        file_object.write(',')
	        file_object.write(param3)
	        file_object.write(',')
	        file_object.write(param4)
	        file_object.write(',')
	        file_object.write(param5)
	        file_object.write(',')
	        file_object.write(param6)
	        file_object.write(',')
	        file_object.write(param7)
	        file_object.write(',')
	        file_object.write(param8)
	        file_object.write(',')
	        file_object.write(param9)
	        file_object.write(',')
	        file_object.write(param10)
	        file_object.write(',')
	        file_object.write(param11)
	        file_object.write(',')
	        file_object.write(param12)

	def loadfile(self):
	    filename = filedialog.askopenfilename(title = "Select file",filetypes = (("txt files","*.txt"),("all files","*.*")))
	    print(filename)
	    readfile = open(filename, "r")
	    for line in readfile:
	    	Type = line.split(",")

	    	self.Nx.set(Type[0])
	    	self.xleft.set(Type[1])
	    	self.xright.set(Type[2])
	    	self.dx.set(Type[3])
	    	self.xC.set(Type[4])
	    	self.xwide.set(Type[5])
	    	self.A0.set(Type[6])
	    	self.hC.set(Type[7])
	    	self.h0.set(Type[8])
	    	self.hshall.set(Type[9]) 
	    	self.hwide.set(Type[10])
	    	self.tfin.set(Type[11])
	    	self.dt.set(Type[12])

	def save_animate_gif(self):
	    file_name_gif = self.Aniname_gif.get()
	    
	    if len(self.Aniname_gif.get())==0:
	    	messagebox.showwarning("Warning", "Please Give Name File!")
	    else:
	    	animation_gif = self.ani.save(file_name_gif+'.gif', writer='imagemagick', fps=30, bitrate=20)
	    	return animation_gif

	def save_animate_mp4(self):
		Writer = animation.writers['pillow']
		writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
		file_name_mp4 = self.Aniname_mp4.get()

		if len(self.Aniname_mp4.get())==0:
			messagebox.showwarning("Warning", "Please Give Name File!")
		else:
			animation_mp4 = self.ani.save(file_name_mp4+'.png', writer='imagemagick')
			return animation_mp4
	        
	def init_window(self):

		# ===========================================================================#
		# 						PARAMETER input                                      #
		# ===========================================================================#
		self.running = False
		self.ani = None

		self.xleft   = StringVar()
		self.xright  = StringVar()
		# domain
		self.Nx      = StringVar()
		self.dx      = StringVar()
		# self.x  = [0]*(self.Nx+1)
		self.tfin = StringVar()
		self.A0 = StringVar()
		self.xC = StringVar()
		self.xwide = StringVar()
		self.dt = StringVar() 
		self.hC = StringVar()
		self.h0 = StringVar()
		self.hshall = StringVar()
		self.hwide = StringVar()

		self.file_name = StringVar()
		self.ani_name_mp4 = StringVar()
		self.ani_name_gif = StringVar()

		self.master.title("MoSV 1D")
		root.iconbitmap("icon.ico")
		vcmd = (self.register(self.validate_float),'%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')

		# --------------------------------------------#
		#  				PARAMETER INPUT 			  #
		# --------------------------------------------#
		tk.Label(self,text="Simulation of SWE 1D",pady=10, padx=50, font ="Helvetica 24").grid(column=0, row=0,columnspan=2)     
		
		border=LabelFrame(root, text = "", pady=5 , padx=2)
		border.grid(row=0,column=0,padx=15,pady=10,columnspan=2)

		self.frame = LabelFrame(border, text = "", pady=8, padx=8)
		self.frame.grid(row=0, column=0, padx=10,pady=2,sticky=W)

		self.L0 = Label(self.frame, text='Domain', font='Helvetica 12 bold', fg='red')
		self.L0.grid(row=0, column=0,columnspan=2, sticky=W, padx= 10)

		self.L0 = Label(self.frame, text='Nx').grid(row=1, column=0, sticky=W, padx= 22)
		self.E0 = tk.Entry(self.frame, textvariable=self.Nx, validatecommand = vcmd)
		self.E0.grid(row=1, column=1, padx= 10)
		self.E0.insert(0,'100')
		self.S0 = Label(self.frame, text='m/s').grid(row=1, column=2,padx=8)

		self.L0 = Label(self.frame, text='Xmin').grid(row=2, column=0, sticky=W, padx= 22)
		self.E1 = tk.Entry(self.frame, textvariable=self.xleft, validatecommand = vcmd)
		self.E1.grid(row=2, column=1, padx= 10)
		self.E1.insert(0,'0')
		self.S0 = Label(self.frame, text='m').grid(row=2, column=2,padx=8)

		self.L0 = Label(self.frame, text='Xmax').grid(row=3, column=0, sticky=W, padx= 22)
		self.E2 = tk.Entry(self.frame, textvariable=self.xright, validatecommand = vcmd)
		self.E2.grid(row=3, column=1, padx= 10)
		self.E2.insert(0,'1.0')
		self.S0 = Label(self.frame, text='m').grid(row=3, column=2,padx=8)

		self.L3 = Label(self.frame, text='dx').grid(row=4, column=0, sticky=W, padx= 22)
		self.E3 = tk.Entry(self.frame, textvariable=self.dx, validatecommand = vcmd)
		self.E3.grid(row=4, column=1, padx= 10)
		self.E3.insert(0,'0.01')
		self.E3.configure(state='readonly')
		self.S3 = Label(self.frame, text='m').grid(row=4, column=2,padx=8)

		# --------------------------------------------#
		#		 	   Wave Initial condition		  #
		# --------------------------------------------#
		self.frame1 = LabelFrame(border, text="", padx=8, pady=8)
		self.frame1.grid(row=1, column=0, padx=10,pady=2,sticky=W)

		self.L4 = Label(self.frame1, text='Wave', font='Helvetica 12 bold', fg='red')
		self.L4.grid(row=3, column=0, columnspan=2, sticky=W, padx= 10)

		self.L4 = Label(self.frame1, text='xC').grid(row=5, column=0, sticky=W, padx= 25)
		self.E4 = tk.Entry(self.frame1, textvariable=self.xC, validatecommand = vcmd)
		self.E4.grid(row=5, column=1, padx= 10)
		self.E4.insert(0,'0.5')
		self.S4 = Label(self.frame1, text='m').grid(row=5, column=2,padx=8)

		self.L5 = Label(self.frame1, text='xWide').grid(row=6, column=0, sticky=W, padx= 25)
		self.E5 = Entry(self.frame1, textvariable=self.xwide, validatecommand = vcmd)
		self.E5.grid(row=6, column=1, padx= 10)
		self.E5.insert(0,'0.05')
		self.S5 = Label(self.frame1, text='m').grid(row=6, column=2,padx=8)

		self.L6 = Label(self.frame1, text='A0').grid(row=7, column=0, sticky=W, padx= 25)
		self.E6 = Entry(self.frame1, textvariable=self.A0, validatecommand = vcmd)
		self.E6.grid(row=7, column=1, padx= 10)
		self.E6.insert(0,'0.002')
		self.S6 = Label(self.frame1, text='m').grid(row=7, column=2,padx=8)

		# --------------------------------------------#
		# 					Bathymetry	              #
		# --------------------------------------------#
		self.frame2 = LabelFrame(border, padx=8, pady=8)
		self.frame2.grid(row=2,column=0,padx=10,pady=2, sticky=W)

		self.L7 = Label(self.frame2, text='Bathymetry',font='Helvetica 12 bold', fg='red')
		self.L7.grid(row=6, column=0, sticky=W, padx= 10,columnspan=2)

		self.L7 = Label(self.frame2, text='hC').grid(row=8, column=0, sticky=W, padx= 25)
		self.E7 = Entry(self.frame2, textvariable=self.hC, validatecommand = vcmd)
		self.E7.grid(row=8, column=1, padx= 10)
		self.E7.insert(0,'750.0')
		self.S7 = Label(self.frame2, text='m').grid(row=8, column=2,padx=8)

		self.L8 = Label(self.frame2, text='h0').grid(row=9, column=0, sticky=W, padx= 25)
		self.E8 = Entry(self.frame2, textvariable=self.h0, validatecommand = vcmd)
		self.E8.grid(row=9, column=1, padx= 10)
		self.E8.insert(0,'20.0')
		self.S8 = Label(self.frame2, text='m').grid(row=9, column=2,padx=8)

		self.L9 = Label(self.frame2, text='hshall').grid(row=10, column=0, sticky=W, padx= 25)
		self.E9 = Entry(self.frame2, textvariable=self.hshall, validatecommand = vcmd)
		self.E9.grid(row=10, column=1, padx= 10)
		self.E9.insert(0,'1.0')
		self.S9 = Label(self.frame2, text='m').grid(row=10, column=2,padx=8)

		self.L10 = Label(self.frame2, text='hwide').grid(row=11, column=0, sticky=W, padx= 25)
		self.E10 = Entry(self.frame2, textvariable=self.hwide, validatecommand = vcmd)
		self.E10.grid(row=11, column=1, padx= 10)
		self.E10.insert(0,'100.0')
		self.S10 = Label(self.frame2, text='m').grid(row=11, column=2,padx=8)

		# --------------------------------------------#
		# 						Time		  		  #
		# --------------------------------------------#
		self.frame3 = LabelFrame(border, padx=8, pady=8)
		self.frame3.grid(row=3,column=0,padx=10,pady=2, sticky=W)

		self.L15 = Label(self.frame3, text='Time',font='Helvetica 12 bold', fg='red')
		self.L15.grid(row=10, column=0, sticky=W, padx= 10,columnspan=2)

		self.L11 = Label(self.frame3, text='tfin').grid(row=12, column=0, sticky=W, padx= 29)
		self.E11 = Entry(self.frame3, textvariable=self.tfin, validatecommand = vcmd)
		self.E11.grid(row=12, column=1, padx= 15)
		self.E11.insert(0,'2.0')
		self.S11= Label(self.frame3, text='s').grid(row=12, column=2,padx=8)

		self.L12 = Label(self.frame3, text='dt').grid(row=13, column=0, sticky=W, padx= 29)
		self.E12 = Entry(self.frame3, textvariable=self.dt, validatecommand = vcmd)
		self.E12.grid(row=13, column=1, padx= 15)
		self.E12.insert(0,'0.001')
		self.E12.configure(state='readonly')
		self.S12 = Label(self.frame3, text='s').grid(row=13, column=2,padx=8)

		#================================================================================#
		# 									 SAVE & LOAD FILE 							 #
		#================================================================================#
		self.frame5 = LabelFrame(border,text="Save Test Case",padx=4,pady=7)
		self.frame5.grid(row=4, column=0,padx=5, columnspan=2, sticky=W)

		self.Lname = Label(self.frame5, text='Test Case').grid(row=0, column=0, padx= 4)
		self.Ename = Entry(self.frame5, textvariable=self.file_name, width=12)
		self.Ename.grid(row=0, column=1, padx= 4)
		self.Pname = Label(self.frame5, text='.txt').grid(row=0, column=2)
		self.B0 = Button(self.frame5, text='SAVE',width=5,bg='green',command=self.save)
		self.B0.grid(row=0, padx=4,column=3)
		self.B1 = Button(self.frame5, text='LOAD',width=5,bg='yellow',command=self.loadfile)
		self.B1.grid(row=0, padx=4,column=4)

		#================================================================================#
		# 								SAVE ANIMATION                                   #
		#================================================================================#
		self.frame6 = LabelFrame(border,text="Save Animation",padx=4,pady=7)
		self.frame6.grid(row=5, column=0,padx=5, columnspan=2, sticky=W)

		# self.Lname = Label(self.frame6, text='Animation ').grid(row=0, column=0, padx= 4)
		self.Aniname_gif = Entry(self.frame6, width=22)
		self.Aniname_gif.grid(row=0, column=1, padx= 4)
		self.Pname_gif = Label(self.frame6, text='.gif').grid(row=0, column=2)
		self.B0 = Button(self.frame6, text='SAVE',width=11,bg='green',command=self.save_animate_gif)
		self.B0.grid(row=0, padx=4,column=3)

		self.Aniname_mp4 = Entry(self.frame6, width=22)
		self.Aniname_mp4.grid(row=1, column=1, padx= 4)
		self.Pname_mp4 = Label(self.frame6, text=' .mp4').grid(row=1, column=2)
		self.B1 = Button(self.frame6, text='SAVE',width=11,bg='green',command=self.save_animate_mp4)
		self.B1.grid(row=1, padx=4,column=3)


		#================================================================================#
		# 									 BUTTON 									 #
		#================================================================================#

		self.frame7 = LabelFrame(border,padx=8,pady=6)
		self.frame7.grid(row=6, column=0, rowspan=2 ,padx=5,pady=2,sticky=W)

		self.B5 = Button(self.frame7,text='RUN',width=16,bg='green',command=self.on_click)
		self.B5.grid(row=14, padx=6,column=0)
		self.B1 = Button(self.frame7,text='CLEAR',width=16,bg='yellow',command=self.clear_text)
		self.B1.grid(row=14, padx=6,column=1)

		#================================================================================#
		# 								 PLOT & ANIMATION								 #
		#================================================================================#

		self.frame6 = LabelFrame(border,padx=10,pady=10)
		self.frame6.grid(row=0, column=3, rowspan=5, padx=10,pady=5,sticky=W)

root.protocol("WM_DELETE_WINDOW", on_closing)
Window(root)
root.mainloop()
