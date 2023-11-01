//Needs folder called "Images"
import java.util.Scanner;//C:\Users\aweso_suq03xb\Desktop\Codestuff\CodeProjects\Particlesim\ParticleSim.java
import java.lang.Math;
import java.util.ArrayList;
import java.awt.*;
import java.awt.image.BufferedImage;
import javax.swing.*;
import java.io.File;
import java.io.*;
import javax.imageio.ImageIO;
public class ParticleSim {
  public static void main(String[] args) {
    Scanner sc = new Scanner(System.in);
    System.out.print("Output Video Length: ");
    double time = sc.nextDouble();
    System.out.print("\nSpeeding up(>1)/Slowing down(<1)/Real Time(1): ");//You can speed it up if necessary by this factor.
    double timerate=sc.nextDouble();
    System.out.print("\nFPS: ");//Controls FPS
    double fpss=sc.nextDouble();
    System.out.print("\nSteps per second (in simulation time): ");
    double step=sc.nextDouble();
    System.out.print("\nParticle Amount: ");
    int num=sc.nextInt();
    System.out.print("\nInitial max vertical speed: ");
    double velx=sc.nextDouble();
    System.out.print("\nInitial max horizontal speed: ");
    double vely=sc.nextDouble();
    System.out.print("\nScreen Width: ");
    double width=sc.nextDouble();
    System.out.print("\nScreen Height: ");
    double height = sc.nextDouble();
    System.out.print("\nParticle Size: ");
    double size=sc.nextDouble();
    System.out.print("\n0 for infinite boundaries w/ PBC, 1 for there to be actual boundaries: ");
    int type = sc.nextInt();
    System.out.println("\nInitial Temperature: ");
    double tempstar=sc.nextDouble();
    System.out.print("\nVideo destination file: ");
    sc.nextLine();
    String path = sc.nextLine();
    double vidlength=time;
    time=time*timerate;
    JFrame framer = new JFrame();
    framer.setSize((int) width,(int) height);
    Particle[] particles=new Particle[num];
    Particle.screenwidth=width;
    Particle.screenheight=height;
    Particle.type=type;
    int thisworks = (int)(Math.log10(step*time*timerate+1)+1);;
    BufferedImage screenshot = new BufferedImage((int) width,(int) height,BufferedImage.TYPE_INT_ARGB);
    for(int i=0; i<num; i++){
      particles[i]=new Particle(size+Math.random()*(width-2*size),size+Math.random()*(height-2*size),Math.random()*(velx)-(velx/2),Math.random()*(vely)-(vely/2),1/step,size,1); //initializes particle array
    }
    Particle.preset();//Sets up particles so none of them are intersecting.
    double k33 = ((Particle.total_kinetic_energy()/Particle.Ledger.size()))/Math.pow(10,23)/1.3/tempstar; //Defines epsilon in Lennard Jones potential using temperature, takes into account boltzmann constant.
    double k34=Particle.total_energy();//Finds the total initial energy of the particles.
    for(int i=0;i<num;i++){
      particles[i].epsilon=k33;//Sets up epsilon in the Lennard-Jones potential through the given initial temperature.
    }
    screenshot=redraw((int) width, (int) height, particles);
    String curdir = System.getProperty("user.dir");
    for(int k=0;k<step*time;k++){//timeframe
      for(int i=0; i<num; i++){//takestep
        particles[i].takestep();
      }
      for(int i=0; i<num; i++){//implement takestep
          particles[i].implement_takestep();
      }
    //System.out.println(Particle.total_energy()+","+Particle.total_momentum()+"\n"+particles[0].pos[0]+", "+particles[0].pos[1]+"\n"+Particle.total_kinetic_energy()+", "+Particle.calc_total_pot_energy()+"\n\n");
    screenshot=redraw((int) width, (int) height, particles);
    if(!(((int)(k*fpss/(step*timerate))-((int)((k-1)*fpss/(step*timerate)))==0))){
      downloadImage(screenshot,(int)(k*fpss/(step*timerate)),thisworks);//Downloads frames in folder called "Image"
    }
  }
    //System.out.println("\n"+(Particle.total_energy()-k34)+"\n"+Particle.total_energy());
    ProcessBuilder asdf = new ProcessBuilder("C:\\Users\\aweso_suq03xb\\Downloads\\ffmpeg-6.0-essentials_build\\ffmpeg-6.0-essentials_build\\bin\\ffmpeg.exe","-y","-r",fpss+"","-f","image2","-s",((int) width)+"x"+((int) height),"-i",curdir+"\\images\\tempoutput%0"+thisworks+"d.png","-vcodec","libx264","-crf","5","-pix_fmt","yuv420p", curdir+"\\"+path);//Compiles image into video.
    System.out.print("Complete.");
    try{
    Process b = asdf.start();
    }
    catch (Exception e){
    e.printStackTrace();
    }
  }
  public static BufferedImage redraw(int width, int height, Particle[] particles){
    BufferedImage a = new BufferedImage(width,height, BufferedImage.TYPE_INT_ARGB);
    Graphics2D gp = a.createGraphics();
    gp.setStroke(new BasicStroke(1));
    for(int i=0;i<particles.length;i++){
      gp.fillOval((int) (particles[i].pos[0]-particles[i].size),(int) (particles[i].pos[1]-particles[i].size),2*((int) particles[i].size),2*((int) particles[i].size));
    }
    gp.dispose();
    return a;
  }//Sets up particles
  public static void downloadImage(BufferedImage a, int n,int j){
    String curdir = System.getProperty("user.dir");
    int s;
    if(n==0){
      s=1;
    }
    else{
      s = (int) (Math.log10(n)+1);
    }
    String aa=new String(new char[(j-s)]).replace('\0','0');
    File newfile = new File(curdir+"\\images\\tempoutput"+aa+n+".png");
    try{
      newfile.createNewFile();
    }
    catch(Exception e){
      System.err.println(e);
    }
    try{
      ImageIO.write(a,"png",newfile);
    }
    catch(Exception e){
      System.err.println(e);
    }
  }//Downloads each image
}
class Particle{
  double[] pos=new double[2];
  double[] vel=new double[2];
  double energy;
  double size;
  double epsilon;
  double listpos;
  double density;
  double timestep;
  double mass;
  static double screenwidth;
  static double screenheight;
  static int type;
  double[] newpos=new double[2];
  double[] newvel=new double[2];
  static ArrayList<Particle> Ledger = new ArrayList<Particle>();
  public Particle(double posx,double posy, double velx, double vely, double timestep, double size, double density){
    this.pos[0]=((posx%screenwidth)+screenwidth)%screenwidth;
    this.pos[1]=((posy%screenheight)+screenheight)%screenheight;
    this.vel[0]=velx;
    this.vel[1]=vely;
    this.listpos=Ledger.size();
    this.timestep=timestep;
    this.size=size;
    this.density= density;
    this.mass=Math.pow(size,2)*Math.PI*density;
    Ledger.add(this);
  }
  public static void preset(){//Sets up particles so none are overlapping
    int d=1;
    double xdiff;
    double ydiff;
    double total;
    double xdis;
    double ydis;
    while (!(d==0)){
      d=0;
      for(int i=0;i<Particle.Ledger.size();i++){
      Particle as = Ledger.get(i);
      for(int j =0;i<(as.listpos);j++){
        Particle s = Ledger.get(j);
        if(Math.sqrt(Math.pow((as.newpos[0]-Ledger.get(i).newpos[0]),2)+Math.pow((as.newpos[1]-Ledger.get(i).newpos[1]),2))<(as.size+Ledger.get(i).size)){
          d=d+1;
          as.newvel[0]=(as.mass-s.mass)*as.newvel[0]/(as.mass+s.mass)+(2*s.mass)*s.newvel[0]/(s.mass+as.mass);
          as.newvel[1]=(as.mass-s.mass)*as.newvel[1]/(as.mass+s.mass)+(2*s.mass)*s.newvel[1]/(s.mass+as.mass);
          s.newvel[0]=(2*as.mass)*as.newvel[0]/(as.mass+s.mass)+(s.mass-as.mass)*s.newvel[0]/(s.mass+as.mass);
          s.newvel[1]=(2*as.mass)*as.newvel[1]/(as.mass+s.mass)+(s.mass-as.mass)*s.newvel[1]/(s.mass+as.mass);
          xdiff=as.newpos[0]-s.newpos[0];
          ydiff=as.newpos[1]-s.newpos[1];
          total=Math.sqrt(Math.pow(xdiff,2)+Math.pow(ydiff,2));
          xdis=((s.size+as.size)*xdiff/(total))-xdiff;
          ydis=((as.size+s.size)*ydiff/total)-ydiff;
          as.newpos[0]=as.newpos[0]+(xdis)*s.mass/(as.mass+s.mass);
          as.newpos[1]=as.newpos[1]+(ydis)*s.mass/(as.mass+s.mass);
          s.newpos[0]=s.newpos[0]+(-xdis)*as.mass/(as.mass+s.mass);
          s.newpos[1]=s.newpos[1]+(-ydis)*as.mass/(as.mass+s.mass);
        }
      }
    }
    }
  }
  public void implement_takestep(){//Actually takes the next step
    this.pos[0]=this.newpos[0];
    this.pos[1]=this.newpos[1];
    this.vel[0]=this.newvel[0];
    this.vel[1]=this.newvel[1];
  }
  public void takestep(){//Calculations to take next step
    double[] es = new double[(int) this.listpos];
    if(type==0){
      double xdist = timestep*this.vel[0];
      double ydist = timestep*this.vel[1];
      this.newpos[0]=(((this.pos[0]+xdist)%this.screenwidth)+this.screenwidth)%this.screenwidth;
      this.newpos[1]=(((this.pos[1]+ydist)%this.screenheight)+this.screenheight)%this.screenheight;
      double[] temp_var=calc_acc();
      this.newvel[0]=this.vel[0]+temp_var[0]*timestep/(Math.sqrt(Math.pow(xdist,2)+Math.pow(ydist,2)));
      this.newvel[1]=this.vel[1]+temp_var[1]*timestep/(Math.sqrt(Math.pow(xdist,2)+Math.pow(ydist,2)));
    }
    else if(type==1){
      double xdist = timestep*this.vel[0];
      double ydist = timestep*this.vel[1];
      this.newpos[0]=(this.pos[0]+xdist);
      this.newpos[1]=(this.pos[1]+ydist);
      double[] temp_var= calc_acc();
      this.newvel[0]=this.vel[0]+temp_var[0]*timestep/(Math.sqrt(Math.pow(xdist,2)+Math.pow(ydist,2)));
      this.newvel[1]=this.vel[1]+temp_var[1]*timestep/(Math.sqrt(Math.pow(xdist,2)+Math.pow(ydist,2)));
      if(this.newpos[0]<0){
        this.newvel[0]=-this.newvel[0];
        this.newpos[0]=-this.newpos[0];
      }
      if(this.newpos[0]>this.screenwidth){
        this.newvel[0]=-this.newvel[0];
        this.newpos[0]=2*this.screenwidth-this.newpos[0];
      }
      if(this.newpos[1]<0){
        this.newvel[1]=-this.newvel[1];
        this.newpos[1]=-this.newpos[1];
      }
      if(this.newpos[1]>this.screenheight){
        this.newvel[1]=-this.newvel[1];
        this.newpos[1]=2*this.screenheight-this.newpos[1];
      }
    }
    for(int i =0;i<(this.listpos);i++){
      Particle s = Ledger.get(i);
      if(Math.sqrt(Math.pow((this.newpos[0]-Ledger.get(i).newpos[0]),2)+Math.pow((this.newpos[1]-Ledger.get(i).newpos[1]),2))<(this.size+Ledger.get(i).size)){
        this.newvel[0]=(this.mass-s.mass)*this.newvel[0]/(this.mass+s.mass)+(2*s.mass)*s.newvel[0]/(s.mass+this.mass);
        this.newvel[1]=(this.mass-s.mass)*this.newvel[1]/(this.mass+s.mass)+(2*s.mass)*s.newvel[1]/(s.mass+this.mass);
        s.newvel[0]=(2*this.mass)*this.newvel[0]/(this.mass+s.mass)+(s.mass-this.mass)*s.newvel[0]/(s.mass+this.mass);
        s.newvel[1]=(2*this.mass)*this.newvel[1]/(this.mass+s.mass)+(s.mass-this.mass)*s.newvel[1]/(s.mass+this.mass);
        double xdiff=this.newpos[0]-s.newpos[0];
        double ydiff=this.newpos[1]-s.newpos[1];
        double total=Math.sqrt(Math.pow(xdiff,2)+Math.pow(ydiff,2));
        double xdis=((s.size+this.size)*xdiff/(total))-xdiff;
        double ydis=((this.size+s.size)*ydiff/total)-ydiff;
        this.newpos[0]=this.newpos[0]+(xdis)*s.mass/(this.mass+s.mass);
        this.newpos[1]=this.newpos[1]+(ydis)*s.mass/(this.mass+s.mass);
        s.newpos[0]=s.newpos[0]+(-xdis)*this.mass/(this.mass+s.mass);
        s.newpos[1]=s.newpos[1]+(-ydis)*this.mass/(this.mass+s.mass);
      }
    if(type==0){//Added in last, if problem check this part of the code.
      this.newpos[0]=((this.newpos[0])%this.screenwidth+this.screenwidth)%this.screenwidth;
      s.newpos[0]=((s.newpos[0])%screenwidth+this.screenwidth)%this.screenwidth;
      this.newpos[1]=((this.newpos[1]%screenheight)+this.screenheight)%this.screenheight;
      s.newpos[1]=((s.newpos[1]%screenheight)+this.screenheight)%this.screenheight;
    }
  }
}
  public double[] calc_acc(){//Uses calc_acc_part to calculate acceleration for each particle
    double[] force= new double[2];
    double[] temp = new double[2];
    for(int i=0;i<Ledger.size();i++){
      if(!(i==this.listpos)){
        temp=calc_acc_part(this,Ledger.get(i));
        force[0]=force[0]+temp[0];
        force[1]=force[1]+temp[1];
      }
    }
    temp[0]=force[0]/(this.mass);
    temp[1]=force[1]/this.mass;
    return temp;
  }
  public double[] calc_acc_part(Particle Particle1, Particle Particle2){//Calculates acceleration for each direction for each particle based on lennard-jones
    double dist = (Math.sqrt((Math.pow((Particle2.pos[0]-Particle1.pos[0]),2)+Math.pow((Particle2.pos[1]-Particle1.pos[1]),2)))/*-Particle1.size-Particle2.size*/);
    double squared=Math.pow(Particle1.size,2)/((Math.pow(dist,2)));
    if(squared<.0001){
      return new double[] {0,0};
    }
    double sixth=Math.pow(squared,3);
    double twelfth=Math.pow(sixth,2);
    double f=24*Particle1.epsilon*(2*twelfth-sixth)*(-1)/dist;
    double[] s = new double[2];
    s[0]=f*(Particle2.pos[0]-Particle1.pos[0])/dist;
    s[1]=f*(Particle2.pos[1]-Particle1.pos[1])/dist;
    return s;
  }

  //Below this point are methods I used to make sure that the system was working by calculating some of its properties.
  public static double calc_Lennard_Jones(Particle Particle1, Particle Particle2){//Calculates lennard-jones potential for each particle
    double squared=Math.pow(Particle1.size,2)/((Math.pow((Particle2.pos[0]-Particle1.pos[0]),2)+Math.pow((Particle2.pos[1]-Particle1.pos[1]),2)));
    if(squared<.0001){
      return 0;
    }
    double sixth=Math.pow(squared,3);
    double twelfth=Math.pow(sixth,2);
    return 4*Particle1.epsilon*(twelfth-sixth);
  }
  public static double calc_pot_energy(Particle Particle1){//Calculates potential energy
    double  energy1=0;
    for(int i=0;i<Ledger.size();i++){
      if(!(i==Particle1.listpos)){
        energy1=energy1+calc_Lennard_Jones(Particle1,Ledger.get(i));
      }
    }
    return energy1;
  }
  public static double calc_total_pot_energy(){//Calculates total potential energy
    double energy1=0;
    for(int i=0;i<Ledger.size();i++){
      energy1=energy1+calc_pot_energy(Ledger.get(i));
    }
    return energy1;
  }
  public static double kinetic_energy(Particle Particle1){//Calculates kinetic energy of one particle
    return Particle1.mass*(Math.pow(Particle1.vel[1],2)+Math.pow(Particle1.vel[0],2))/2;
  }
  public static double total_kinetic_energy(){//Calculates total kinetic energy
    double energy1=0;
    for(int i=0;i<Ledger.size();i++){
      energy1=energy1+kinetic_energy(Ledger.get(i));
    }
    return energy1;
  }
  public static double total_energy(){//Checks the total energy of the system
    return calc_total_pot_energy()+total_kinetic_energy();
  }
  public static double momentum(Particle Particle1){//Calculates momentum of one particle
    return Particle1.mass*Math.sqrt(Math.pow(Particle1.vel[0],2)+Math.pow(Particle1.vel[1],2));
  }
  public static double total_momentum(){//Calculates total momentum.
    double energy1=0;
    for(int i=0;i<Ledger.size();i++){
      energy1=energy1+momentum(Ledger.get(i));
    }
    return energy1;
  }
}
