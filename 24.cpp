/*
This script uses the second set of equations from here: https://en.wikipedia.org/wiki/Leapfrog_integration#Algorithm
and treats the system as an n-body problem
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <chrono>
#include <fstream>
#include <sstream>

class Particle {
public:
  double x, y, vx, vy, ax, ay, PE, n, tp;
  //double v = sqrt((vx*vx) + (vy*vy)); //this just initialises v but doesn't update it when vx and vy change
  Particle(double nx, double ny, double nvx, double nvy): x(nx), y(ny), vx(nvx), vy(nvy) {}
  double v() {
    return sqrt((vx*vx) + (vy*vy));
  }
  double a() {
    return sqrt((ax*ax) + (ay*ay));
  }
  double KE() { //kinetic energy per unit mass
    return (1.0/2.0)*(v()*v());
  }
  double r() { //distance from origin
    return sqrt((x*x) + (y*y));
  }
  double phi() {
    return atan2(y,x);
  }
  double phic(double phib) {
    return phi() - phib;
  }
  double xc(double phib) {
    return r()*cos(phic(phib));
  }
  double yc(double phib) {
    return r()*sin(phic(phib));
  }
  double vr() {
    return (vy*(y/r())) + (vx*(x/r()));
  }
  double vphi() {
    return (vy*(x/r())) - (vx*(y/r()));
  }
  double phidot() {
    return vphi()/r();
  }
  double vxc(double phib, double omega) {
    return vr()*cos(phic(phib)) - r()*phidot()*sin(phic(phib));
  }
  double vyc(double phib, double omega) {
    return vr()*sin(phic(phib)) + r()*phidot()*cos(phic(phib));
  }
  double L() { //angular momentum per unit mass
    return r()*vphi();
  }
  double E() {
    return KE() + PE;
  }
  double Ej(double omega) {
    return E() - omega*L();
  }
  void print(std::ofstream &ostr) {//, double phib, double omega) { //print all its parameters comma-separated
    //ostr << "," << n << "," << m << "," << x << "," << y << "," << r() << "," << phi() << "," << xc(phib) << "," << yc(phib) << "," << vx << "," << vy << "," << vr() << "," << vphi() << "," << v() << "," << L() << "," << ax << "," << ay << "," << a() << "," << KE() << "," << PE << "," << E() << "," << Ej(omega);
    ostr << "," << n << "," << x << "," << y << "," << vx << "," << vy << "," << ax << "," << ay << "," << PE;
  }
};

void getfromfile(std::vector<Particle> &ps, double &t_start, double &t_end, std::vector<std::string> &psn, bool sos, double &Rcr, double &omega, double &phib0) {
    std::ifstream file("init.txt");
    std::string headers, header, values, value;
    std::getline(file, headers);
    //std::cout << headers << "\n";
    std::stringstream headersstream(headers);

    std::getline(file, values);
    //std::cout << values << "\n";
    std::stringstream valuesstream(values);

    file.close();

    while (std::getline(headersstream, header, ',')) {
        std::size_t pos = header.find('_');
        if (pos != std::string::npos) {
            //std::cout << pos << "\n";
            std::string num = header.substr(pos+1,header.length()-pos);
            //std::cout << strnum << "\n";
            //int num = std::stoi(strnum);
            //std::cout << num << "\n";
            if (!psn.empty()) {
                if (num != psn.back()) {
                    psn.push_back(num);
                }
                //std::cout << psn.back() << "\n";
            } else {
                psn.push_back(num);
            }
        }
    }

    //std::cout << psn.size() << "\n";

    /*
    for (int i=0; i<psn.size(); i++) {
        std::cout << psn[i] << "\n";
    }
    */

    std::vector<double> vals;
    while (std::getline(valuesstream, value, ',')) {
        double val = std::stod(value);
        vals.push_back(val);
    }
    t_start = vals[0];
    t_end = vals[1];
    int i_start;
    if (sos) {
        phib0 = vals[2];
        omega = vals[3];
        Rcr = vals[4];
        i_start = 5;
    } else {
        i_start = 2;
    }
    for (int i = i_start; i < vals.size(); i += 4) {
        //std::cout << vals[i] << vals[i+1] << vals[i+2] << vals[i+3] << "\n";
        ps.push_back(Particle(vals[i],vals[i+1],vals[i+2],vals[i+3]));
    }
}

void stuff(double &t, double t1, double &A, double Af, double &phib, double vc2, double &Rcr, double &omega, double Rcr0, double t2, double vcr, double a, double phib1, double phib2, double eta, double omega0, bool sos, double phib0, double t_start) {

  if (sos) {
      A = Af;
      phib = phib0 + omega*(t - t_start);
  } else {
      if (t <= t1) {
        double xi = 2.0*(t/t1)-1;
        A = Af*(0.1875*pow(xi,5) - 0.625*pow(xi,3) + 0.9375*xi + 0.5);
        omega = sqrt(vc2)/Rcr0;
        phib = omega*t;
      } else if (t > t2) {
        A = Af;
        Rcr = Rcr0 + 0.5*vcr*(t2-t1) + vcr*(t-t2);
        omega = sqrt(vc2)/(Rcr);
        phib = (1.0/eta)*log(Rcr/(Rcr0 + 0.5*vcr*(t2-t1))) + phib2 + phib1;
      } else if (t > t1) {
        A = Af;
        Rcr = Rcr0 + 0.5*vcr*((t-t1)*(t-t1)/(t2-t1));
        omega = sqrt(vc2)/(Rcr);
        phib = (2.0*(t2-t1)/(eta*a))*atan((t-t1)/a) + phib1;
      }

      /*
      Rcr = Rcr0;
      A = Af;
      omega = omega0;
      phib = omega0*t;
      */
  }
}

void timestep(double tp, double k, double a, double &n, double p, double r, double PE) {
  double a1 = log2(1.0/(k/pow(a,p)));
  double tn;
  if (a1 < 0) {
    tn = 0;
  } else {
    tn = ceil(a1);
  }
  //std::cout << "\n" << tn << "," << a1 << "," << PE << "," << a;

  if (n > tn) {
      for (double d = n; d >= tn; d -= 1) {
          double dtf = 1.0/pow(2,d);
          double mot = tp/dtf;
          if (floor(mot) == mot) {
            n = d;
          }
      }
  } else {
    n = tn;
  }

}

void acc(std::vector<Particle> &ps, int i, double G, double vc2, double A, int M, double Rcr, double b, double phib, double kpc, double Rb) { //calculate net acceleration of particle i
  double R = ps[i].r(), phi = ps[i].phi();
  ps[i].PE = vc2*log(R/kpc) - (A*vc2/M)*pow((R/Rcr),2)*pow((b+1)/(b+R/Rcr),5)*cos(M*(phi-phib));
  double abarr = (1.0/Rcr)*(A*vc2/M)*(R/Rcr)*pow((b+1)/(b+(R/Rcr)),5)*(2.0-(5.0/((b*Rcr/R)+1.0)))*cos(M*(phi-phib));
  double abarphi = -(A*vc2/R)*pow((R/Rcr),2)*pow((b+1)/(b+(R/Rcr)),5)*sin(M*(phi-phib));
  //ps[i].PE = vc2*log(R/kpc) - (A*vc2/M)*pow((R/Rcr),2)*pow((Rb+Rcr)/(Rb+R),5)*cos(M*(phi-phib));
  //double abarr = (1.0/Rcr)*(A*vc2/M)*(R/Rcr)*pow((Rb+Rcr)/(Rb+R),5)*(2.0-(5.0/((Rb/R)+1.0)))*cos(M*(phi-phib));
  //double abarphi = -(A*vc2/R)*pow((R/Rcr),2)*pow((Rb+Rcr)/(Rb+R),5)*sin(M*(phi-phib));
  ps[i].ax = -vc2*ps[i].x/(R*R) + abarr*(ps[i].x/R) - abarphi*(ps[i].y/R);
  ps[i].ay = -vc2*ps[i].y/(R*R) + abarr*(ps[i].y/R) + abarphi*(ps[i].x/R);
}

int main() {
  //std::ostringstream ostr;
  std::ofstream ostr;
  ostr.open("data10.txt");
  auto start = std::chrono::high_resolution_clock::now();

  double Myr = 3.154e13, km = 1.0e3, kpc = 3.086e19;
  double t_start = 0.0, dtmax = 1.0e-1*Myr/6.0, p = 1.0, t_end = Myr*20e3;
  double k = pow(7.18e-9, p);
  double G = 6.67e-11, vc2 = pow(235*km,2), A = 0.0, Af = 0.02, eta = 0.003, omega = 80*km/kpc, b = 0.75, t1 = Myr*2e3, phib = 0;
  double t2 = t1 + (Myr*1e3), Rcr0 = sqrt(vc2)/omega;
  double Rcr = Rcr0, Rb = b*Rcr0;
  double vcr = sqrt(vc2)*eta;
  int M = 2;
  double ax_old, ay_old;
  double omega0 = omega;
  double tf = 0;
  double phib0;

  double x0 = -1.0/(pow(2.0,2.0/3.0) - 1), x1 = 1.0/(2.0 - pow(2.0,1.0/3.0));
  double c1 = 0.5*x1, c2 = 0.5*(x0+x1), c3 = 0.5*(x0+x1), c4 = 0.5*x1, d1 = x1, d2 = x0, d3 = x1;
  std::vector<double> coefficients = {c1,c2,c3,d1,d2,d3};
  //std::cout << (coefficients[0] == c1) << "," << (coefficients[1] == c2) << "," << (coefficients[2] == c3) << "," << (coefficients[3] == d1) << "," << (coefficients[4] == d2) << "," << (coefficients[5] == d3);

  std::vector<Particle> ps;

  /*
  std::default_random_engine eng;
  eng.seed(start.time_since_epoch().count());
  std::uniform_real_distribution<double> disphi(0.0,2.0*M_PI);
  std::uniform_real_distribution<double> disr(b*Rcr,2.0*b*Rcr);
  std::uniform_real_distribution<double> dism(1.0e1,1.0e2);
  for (int i=1; i <= 60; i++) {
    double r0 = 1.9*b*Rcr;
    double n = i

    if (i==20) {
      r0 = 0.33*Rcr;
    }

    //double phi0 = 0;// disphi(eng);
    //double m = dism(eng);
    ps.push_back(Particle(1.0, r0, 0, 0, sqrt(vc2), n));
  }
  */



  double r_start = 0.6*Rcr;
  double r_end = 3.0*Rcr;
  double dr = 0.05*Rcr;
  int nphis = 24;
  double dphi0 = 0;
  double phi0 = 0;
  for (double r = r_start; r <= r_end; r += dr) {
    for (int i=0; i<nphis; i++) {
        double phi = phi0 + i*2*M_PI/nphis;
        ps.push_back(Particle(r*cos(phi), r*sin(phi), -sqrt(vc2)*sin(phi), sqrt(vc2)*cos(phi)));
    }
    phi0 += dphi0;
  }



/*
  ps = {
    Particle(1.0*Rcr,0,0,sqrt(vc2)),
    Particle(1.05*Rcr,0,0,sqrt(vc2)),
    Particle(1.1*Rcr,0,0,sqrt(vc2)),
    Particle(1.15*Rcr,0,0,sqrt(vc2)),
    Particle(1.2*Rcr,0,0,sqrt(vc2)),
    Particle(1.25*Rcr,0,0,sqrt(vc2)),
    Particle(1.3*Rcr,0,0,sqrt(vc2)),
    Particle(1.35*Rcr,0,0,sqrt(vc2)),
    Particle(1.4*Rcr,0,0,sqrt(vc2)),
    Particle(1.45*Rcr,0,0,sqrt(vc2)),
    Particle(1.5*Rcr,0,0,sqrt(vc2)),
    Particle(1.55*Rcr,0,0,sqrt(vc2)),
    Particle(1.6*Rcr,0,0,sqrt(vc2)),
    Particle(1.65*Rcr,0,0,sqrt(vc2)),
    Particle(1.7*Rcr,0,0,sqrt(vc2)),
    Particle(1.75*Rcr,0,0,sqrt(vc2)),
    Particle(1.8*Rcr,0,0,sqrt(vc2)),
    Particle(1.85*Rcr,0,0,sqrt(vc2)),
    Particle(1.9*Rcr,0,0,sqrt(vc2)),
    Particle(1.95*Rcr,0,0,sqrt(vc2)),
  };
*/


  //ps = {Particle(0,0.5*b*Rcr,-1.8*sqrt(vc2),-10000)};

/*
  //double L0 = 2*0.85*Rcr0*sqrt(vc2);
  //double E0 = 4.3*vc2*(0.5 + log(0.85*Rcr0/kpc));

  double r_start = 0.75*Rcr;
  double r_end = r_start;//1.2*Rcr;
  double dr = 0.05*Rcr;
  int nphis = 24;
  double dphi0 = 0;
  double phi0 = 0;
  double vr_start = -0.6;
  double vr_end = 0.6;
  double dvr = 0.05;
  for (double r = r_start; r <= r_end; r += dr) {
    for (int i=0; i<nphis; i++) {
        double phi = phi0 + i*2*M_PI/nphis;
        for (double fvr = vr_start; fvr <= vr_end; fvr += dvr) {
            double vr = fvr*sqrt(vc2);
            ps.push_back(Particle(r*cos(phi), r*sin(phi), vr*cos(phi) -sqrt(vc2)*sin(phi), vr*sin(phi) + sqrt(vc2)*cos(phi)));
        }
        //double phi = phi0 + i*2*M_PI/nphis;
        //double vphi0 = L0/r;
        //double vr0 = sqrt((2*E0) - (vphi0*vphi0) -2*vc2*log(r/kpc));
        //ps.push_back(Particle(r*cos(phi), r*sin(phi), vr0*cos(phi) -vphi0*sin(phi), vr0*sin(phi) + vphi0*cos(phi)));
        //vr0 = -vr0; //taking the negative square root as well
        //ps.push_back(Particle(r*cos(phi), r*sin(phi), vr0*cos(phi) -vphi0*sin(phi), vr0*sin(phi) + vphi0*cos(phi)));
    }
    phi0 += dphi0;
  }
*/

/*
  double r_start = 0.75*Rcr;
  double r_end = r_start;//1.2*Rcr;
  double dr = 0.05*Rcr;
  int nphis = 24;
  double dphi0 = 0;
  double phi0 = 0;
  double exv_start = 0;
  double exv_end = 0.6;
  double dexv = 0.05;
  for (double r = r_start; r <= r_end; r += dr) {
    for (int i=0; i<nphis; i++) {
        double phi = phi0 + i*2*M_PI/nphis;
        for (double fexv = exv_start; fexv <= exv_end; fexv += dexv) {
            double vphi = sqrt(vc2*(1+(fexv*fexv)));
            ps.push_back(Particle(r*cos(phi), r*sin(phi), -vphi*sin(phi), vphi*cos(phi)));
        }
        //double phi = phi0 + i*2*M_PI/nphis;
        //double vphi0 = L0/r;
        //double vr0 = sqrt((2*E0) - (vphi0*vphi0) -2*vc2*log(r/kpc));
        //ps.push_back(Particle(r*cos(phi), r*sin(phi), vr0*cos(phi) -vphi0*sin(phi), vr0*sin(phi) + vphi0*cos(phi)));
        //vr0 = -vr0; //taking the negative square root as well
        //ps.push_back(Particle(r*cos(phi), r*sin(phi), vr0*cos(phi) -vphi0*sin(phi), vr0*sin(phi) + vphi0*cos(phi)));
    }
    phi0 += dphi0;
  }
*/

  std::vector<std::string> psn;

  bool sos = false; //true if it's for particles with the same Ej for surfaces of section, otherwise false
  //getfromfile(ps,t_start,t_end,psn,sos,Rcr,omega,phib0);

  if (psn.empty()) {
      for (int i=0; i < ps.size(); i++) {
        psn.push_back(std::to_string(i+1));
      }
  }


  double t = t_start;

  double phib1 = (sqrt(vc2)/Rcr0)*t1;
  double a = sqrt(2.0*Rcr0*(t2-t1)/vcr);
  double phib2 = (2.0*(t2-t1)/(eta*a))*atan((t2-t1)/a);
  stuff(t,t1,A,Af,phib,vc2,Rcr,omega,Rcr0,t2,vcr,a,phib1,phib2,eta,omega0,sos,phib0,t_start);
  //Particle baro(0,b*Rcr*cos(phib),b*Rcr*sin(phib),0,0);

  //calculate initial acceleration for all particles:
  for (int i=0; i < ps.size(); i++) {
    acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc,Rb);
    ps[i].n = 0; //initialising as 0, else it could be initialised as a small but non-zero number depending on compiler
    ps[i].tp = 0;
    timestep(ps[i].tp,k,ps[i].a(),ps[i].n,p,ps[i].r(),ps[i].PE); //just calculating the timestep here in order to output
  }

  ostr << "t";
  for (int i=0; i < psn.size(); i++) {
    ostr << ",n_" << psn[i] << ",x_" << psn[i] << ",y_" << psn[i] << ",vx_" << psn[i] << ",vy_" << psn[i] << ",ax_" << psn[i] << ",ay_" << psn[i] << ",PE_" << psn[i];
  }
  ostr << ",phib,omega,Rcr\n" << t;
  for (int i=0; i < ps.size(); i++) {
    ps[i].print(ostr);//, phib, omega);
  }
  ostr << "," << phib << "," << omega << "," << Rcr;
  ostr << "\n";

  while (t < t_end) {
    //double f = 1.0/pow(2,maxn);
    t += dtmax;
    tf += 1;
    /*
    bool cond = (floor(outint) == outint);
    if (cond) {
      ostr << t;
    }
    */
    for (int i=0; i < ps.size(); i++) {
      ps[i].tp = 0;
      while (ps[i].tp < 1.0) {
        double dtf = 1.0/(pow(2,ps[i].n));
        /*
        double dt;
        if (dtf > 1 - ps[i].tp) {
            dt = 1 - ps[i].tp;
        } else {
            dt = dtf*dtmax;
        }
        */
        double dt = dtf*dtmax;
        for (int j=0; j<=2; j++) {
            ps[i].x = ps[i].x + coefficients[j]*ps[i].vx*dt;
            ps[i].y = ps[i].y + coefficients[j]*ps[i].vy*dt;
            acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc,Rb);
            ps[i].vx = ps[i].vx + coefficients[j+3]*ps[i].ax*dt;
            ps[i].vy = ps[i].vy + coefficients[j+3]*ps[i].ay*dt;
        }
        /*
        ps[i].x = ps[i].x + c1*ps[i].vx*dt;
        ps[i].y = ps[i].y + c1*ps[i].vy*dt;
        acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc);
        ps[i].vx = ps[i].vx + d1*ps[i].ax*dt;
        ps[i].vy = ps[i].vy + d1*ps[i].ay*dt;
        ps[i].x = ps[i].x + c2*ps[i].vx*dt;
        ps[i].y = ps[i].y + c2*ps[i].vy*dt;
        acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc);
        ps[i].vx = ps[i].vx + d2*ps[i].ax*dt;
        ps[i].vy = ps[i].vy + d2*ps[i].ay*dt;
        ps[i].x = ps[i].x + c3*ps[i].vx*dt;
        ps[i].y = ps[i].y + c3*ps[i].vy*dt;
        acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc);
        ps[i].vx = ps[i].vx + d3*ps[i].ax*dt;
        ps[i].vy = ps[i].vy + d3*ps[i].ay*dt;
        */
        ps[i].x = ps[i].x + c4*ps[i].vx*dt;
        ps[i].y = ps[i].y + c4*ps[i].vy*dt;
        acc(ps,i,G,vc2,A,M,Rcr,b,phib,kpc,Rb);
        ps[i].tp += dtf;
        timestep(ps[i].tp,k,ps[i].a(),ps[i].n,p,ps[i].r(),ps[i].PE);
      }
      /*
      if (cond) { //checking if the largest integer not greater than outint is equal to outint, in which case it is indeed an integer
        ps[i].print(ostr);//, phib, omega);
      }
      */
    }
    stuff(t,t1,A,Af,phib,vc2,Rcr,omega,Rcr0,t2,vcr,a,phib1,phib2,eta,omega0,sos,phib0,t_start);
    //calculate the new position of the end of the bar:
    /*
    if (cond) {
      ostr << "," << phib << "," << omega << "," << Rcr;
      ostr << "\n";
    }
    */
    double outint = tf/60.0;
    if (floor(outint) == outint) {
        ostr << t;
        for (int i=0; i < ps.size(); i++) {
            ps[i].print(ostr);
        }
        ostr << "," << phib << "," << omega << "," << Rcr << "\n";
    }

  }
  ostr << "\n";
  auto stop = std::chrono::high_resolution_clock::now();
  auto diff = stop - start;
  auto diff_sec = std::chrono::duration_cast<std::chrono::nanoseconds>(diff);
  ostr << "#execution time: " << diff_sec.count() << " nanoseconds";

  /*
  a1 = log2(dtmax/(k/ps[i].a));
  if (a1 < 0) {
    n = 0;
  } else {
    n = ceil(a1);
  }
  dt = dtmax/(pow(2,n));
  for (int i=1; i <= pow(2,n); i++) {
    t += dt;
    ... (calculate acc, integrate)
  }
  */

  /*
  std::string str1 = ostr.str();
  std::ofstream file;
  file.open("data10.txt");
  file << str1;
  file.close();
  */

  ostr.close();
}
