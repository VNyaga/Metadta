/*
Name:	lmetadta.do
Creator:	Victoria N Nyaga
Date: 29th May 2020
Purpose:	Facilitate maintanance of mata functions in lmetadta.mlib.
			These functions are used in metadta.
*/

cd "C:\ado\plus\l"
mata:mata mlib create lmetadta, replace

version 14.0 
mata:
	mata clear
//************************************************************************	
	void koopmancifun(real scalar RR, real scalar y, real rowvector x){
		x1 = x[1]
		n1 = x[2]
		x2 = x[3]
		n2 = x[4]
		level = x[5]
	
		real scalar p1_tilde
		
		ki = invchi2(1, 1 - level)		
		
		p1_tilde = (RR*(n1 + x2) + x1 + n2 - sqrt((RR*(n1 + x2) + x1 + n2)^2 - 4*RR*(n1 + n2)*(x1 + x2)))/(2*(n1 + n2))
		y = (((x1 - n1*p1_tilde)^2)/(n1*p1_tilde*(1 - p1_tilde)))*(1 + (n1*(RR - p1_tilde))/(n2*(1 - p1_tilde))) -  ki
	}

//************************************************************************
//Added on 17 Jan 2019
void koopman_ci(real rowvector v, real scalar alpha) {
	real scalar zstar, x, m, y, n, CIL, CIU, nrat, varhat
	
	zstar = invnormal(1 - alpha/2)
	x = v[1] 
	m = v[2] 
	y = v[3] 
	n = v[4] 
	
	if (x == 0 & y == 0) {
		CIL = 0
		CIU = .
	}
	else {
	    a1 = n * (n * (n + m) * x + m * (n + x) * (zstar^2))
        a2 = -n * (n * m * (y + x) + 2 * (n + m) * y * x + m * (n + y + 2 * x) * (zstar^2))
        a3 = 2 * n * m * y * (y + x) + (n + m) * (y^2) * x + n * m * (y + x) * (zstar^2)
        a4 = -m * (y^2) * (y + x)
        b1 = a2/a1
        b2 = a3/a1
        b3 = a4/a1
        c1 = b2 - (b1^2)/3
        c2 = b3 - b1 * b2/3 + 2 * (b1^3)/27
        ceta = (acos(sqrt(27) * c2/(2 * c1 * sqrt(-c1))))
        t1 = (-2 * sqrt(-c1/3) * cos(pi()/3 - ceta/3))
        t2 = (-2 * sqrt(-c1/3) * cos(pi()/3 + ceta/3))
        t3 = (2 * sqrt(-c1/3) * cos(ceta/3))
        p01 = t1 - b1/3
        p02 = t2 - b1/3
        p03 = t3 - b1/3
        p0sum = p01 + p02 + p03
        p0up = min((p01, p02, p03))
        p0low = p0sum - p0up - max((p01, p02, p03))
        
        rat = (x/m)/(y/n)
        nrat = (x/m)/(y/n)
        varhat = (1/x) - (1/m) + (1/y) - (1/n)
        
        if ((x == 0) & (y != 0)) {
            nrat = ((x + 0.5)/m)/(y/n)
            varhat = (1/(x + 0.5)) - (1/m) + (1/y) - (1/n)
        }
        if ((y == 0) & (x != 0)) {
            nrat = (x/m)/((y + 0.5)/n)
            varhat = (1/x) - (1/m) + (1/(y + 0.5)) - (1/n)
        }
        if ((y == n) & (x == m)) {
            nrat = 1
            varhat = (1/(m - 0.5)) - (1/m) + 1/(n - 0.5) - (1/n)
        }
        La = nrat * exp(-1 * zstar * sqrt(varhat)) * 1/4
        Ha = nrat
	
		if ((x != 0) & (y == 0)) {
		  if (x == m) {
			CIL = (1 - (m - x) * (1 - p0low)/(y + m - (n + m) * p0low))/p0low
			CIU = .
		  }
		  else {
			S = solvenl_init()
			solvenl_init_evaluator(S, &koopmancifun())
			solvenl_init_type(S, "zero")
			solvenl_init_technique(S, "newton")
			solvenl_init_numeq(S, 1)
			solvenl_init_narguments(S, 1)	
			solvenl_init_argument(S, 1, (x, m, y, n, alpha))
			solvenl_init_startingvals(S, La)
			solvenl_init_iter_log(S, "off")
			
			CIL = solvenl_solve(S)
			CIU = .
		  }
		}
	
		if ((x == 0) & (y != n)) {
			S = solvenl_init()
			solvenl_init_evaluator(S, &koopmancifun())
			solvenl_init_type(S, "zero")
			solvenl_init_technique(S, "newton")
			solvenl_init_numeq(S, 1)
			solvenl_init_narguments(S, 1)	
			solvenl_init_argument(S, 1, (x, m, y, n, alpha))
			solvenl_init_startingvals(S, Ha)
			solvenl_init_iter_log(S, "off")
			
			CIU = solvenl_solve(S)
			CIL = 0
		}
	
		if (((x == m) | (y == n)) & (y != 0)) {
		  if ((x == m) & (y == n)) {
			CIL = m/(m + invchi2(1, 1 - alpha))
			CIU = (n + invchi2(1, 1 - alpha))/n	
		  }
		  if ((x == m) & (y != n)) {
			phat1 = x/m
			phat2 = y/n
			phihat = phat2/phat1
			phiu = 1.1 * phihat
			r = 0
			while (r >= -zstar) {
			  a = (m + n) * phiu
			  b = -((x + n) * phiu + y + m)
			  c = x + y
			  p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
			  p2hat = p1hat * phiu
			  q2hat = 1 - p2hat
			  var = (m * n * p2hat)/(n * (phiu - p2hat) + m * q2hat)
			  r = ((y - n * p2hat)/q2hat)/sqrt(var)
			  phiu1 = phiu
			  phiu = 1.0001 * phiu1
			}
			CIU = (1 - (m - x) * (1 - p0up)/(y + m - (n + m) * p0up))/p0up
			CIL = 1/phiu1
		  }
	  
		  if ((y == n) & (x != m)) {
			phat2 = y/n
			phat1 = x/m
			phihat = phat1/phat2
			phil = 0.95 * phihat
			r = 0
			if (x != 0) {
			  while (r <= zstar) {
				a = (n + m) * phil
				b = -((y + m) * phil + x + n)
				c = y + x
				p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
				p2hat = p1hat * phil
				q2hat = 1 - p2hat
				var = (n * m * p2hat)/(m * (phil - p2hat) + n * q2hat)
				r = ((x - m * p2hat)/q2hat)/sqrt(var)
				CIL = phil
				phil = CIL/1.0001
			  }
			}
			
			phiu = 1.1 * phihat
			if (x == 0) {
			  CIL = 0
			  if (n < 100) {
				phiu =  0.01
			  }
			  else {
				phiu = 0.001
			  }
			}
			
			r = 0
			while (r >= -zstar) {
			  a = (n + m) * phiu
			  b = -((y + m) * phiu + x + n)
			  c = y + x
			  p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
			  p2hat = p1hat * phiu
			  q2hat = 1 - p2hat
			  var = (n * m * p2hat)/(m * (phiu - p2hat) + n * q2hat)
			  r = ((x - m * p2hat)/q2hat)/sqrt(var)
			  phiu1 = phiu
			  phiu = 1.0001 * phiu1
			}
			
			CIU = phiu1
		  }
		}
	
		else if ((y != n) & (x != m) & (x != 0) & (y != 0)) {
		
			inits_l = rat*exp(-1*invnormal(1 - alpha/2)*sqrt(1/x + 1/y - 1/m - 1/n))
			inits_u = rat*exp(invnormal(1 - alpha/2)*sqrt(1/x  + 1/y - 1/n - 1/m))
			
			S = solvenl_init()
			solvenl_init_evaluator(S, &koopmancifun())
			solvenl_init_type(S, "zero")
			solvenl_init_technique(S, "newton")
			solvenl_init_numeq(S, 1)
			solvenl_init_narguments(S, 1)	
			solvenl_init_argument(S, 1, (x, m, y, n, alpha))
			solvenl_init_startingvals(S, inits_l)
			solvenl_init_iter_log(S, "off")
			CIL = solvenl_solve(S)
		  
			S = solvenl_init()
			solvenl_init_evaluator(S, &koopmancifun())
			solvenl_init_type(S, "zero")
			solvenl_init_technique(S, "newton")
			solvenl_init_numeq(S, 1)
			solvenl_init_narguments(S, 1)	
			solvenl_init_argument(S, 1, (x, m, y, n, alpha))
			solvenl_init_startingvals(S, inits_u)
			solvenl_init_iter_log(S, "off")
			CIU = solvenl_solve(S)
		}
	}
	ci = (CIL, CIU)
	
	st_matrix("ci", ci)	
}

	mata mlib add lmetadta *()
	mata mlib index
end
