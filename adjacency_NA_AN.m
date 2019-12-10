function [B B1 B1_r B1_u B2]=adjacency_NA_AN(ncell_N,ncell_A,A, P_B, P_B1, P_B1_r, P_B1_u, P_B2)

%% B for VIP N2A

  for j=1:ncell_N
          random_number = rand(1);
          if(random_number <= P_B)
              lamda_B(j) = 1;
          else
              lamda_B(j)= 0;
          end
  end

B = lamda_B';

%% B1 for GABA A2N


  for p=1:ncell_N
          random_number = rand(1);
          if(random_number <= P_B1)
              lamda_B1(p) = 1;
          else
              lamda_B1(p)= 0;
          end
  end

B1 = lamda_B1';




  for k=1:ncell_N
          random_number = rand(1);
          if(random_number <= P_B1_r)
              lamda_B1_r(k) = 1;
          else
              lamda_B1_r(k)= 0;
          end
  end

B1_r = lamda_B1_r';




  for m=1:ncell_N
          random_number = rand(1);
          if(random_number <= P_B1_u)
              lamda_B1_u(m) = 1;
          else
              lamda_B1_u(m)= 0;
          end
  end

B1_u = lamda_B1_u';


%% B2 for Glu A2N

  for n=1:ncell_N
          random_number = rand(1);
          if(random_number <= P_B2)
              lamda_B2(n) = 1;
          else
              lamda_B2(n)= 0;
          end
  end

B2 = lamda_B2';


