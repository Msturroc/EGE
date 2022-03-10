using Distributions, StatsBase, Statistics # load required packages
# define simple stochastic simulation algorithm for intracellular reactions
function ssa_antibiotic(b1::Float64, q1::Float64, mu::Float64, b2::Float64, q2::Float64, D_1::Float64, K_5::Float64, K_6::Float64, V::Float64, timestep::Float64, d1::Float64, d2::Float64, n::Int64, np::Int64, r::Int64, rp::Int64, NN::Int64, ext_chl::Int64, int_chl::Int64,KK::Int64)
    a = zeros(10)
    a[1] = b1 + q1 * mu
    a[2] = d1 * n
    a[3] = b1 + q1 * mu
    a[4] = d1 * r
    a[5] = (b2 + q2 * mu) * n
    a[6] = (b2 + q2 * mu) * r
    a[7] = D_1 * ext_chl * NN / KK
    a[8] = (K_5 * np / V) / (1 + K_6 / (int_chl / V))
    a[9] = d2 * np
    a[10] = d2 * rp
    asum = sum(a)
    time = 0.0
    while time < timestep
        tau = log(1 / rand()) / asum

        time = time + tau

        drxn = rand() * asum
        psm = 0.0
        j = 0

        while psm < drxn
            j = j + 1
            psm = psm + a[j]
        end

        if time > timestep
            return Int(n), Int(np), Int(r), Int(rp), Int(int_chl), Int(ext_chl)
        end

        if j == 1
            n = n + 1
            asum = asum - a[2] - a[5]
            a[2] = d1 * n
            a[5] = (b2 + q2 * mu) * n
            asum = asum + a[2] + a[5]
        elseif j == 2
            n = n - 1
            asum = asum - a[2] - a[5]
            a[2] = d1 * n
            a[3] = (b2 + q2 * mu) * n
            asum = asum + a[2] + a[5]
        elseif j == 3
            r = r + 1
            asum = asum - a[4] - a[6]
            a[4] = d1 * r
            a[6] = (b2 + q2 * mu) * r
            asum = asum + a[4] + a[6]
        elseif j == 4
            r = r - 1
            asum = asum - a[4] - a[6]
            a[4] = d1 * r
            a[6] = (b2 + q2 * mu) * r
            asum = asum + a[4] + a[6]
        elseif j == 5
            np = np + 1
            asum = asum - a[8] - a[9]
            a[8] = K_5 * ((np / V) / (1 + K_6 / (int_chl / V)))
            a[9] = d2 * np
            asum = asum + a[8] + a[9]
        elseif j == 6
            rp = rp + 1
            asum = asum - a[10]
            a[10] = d2 * rp
            asum = asum + a[10]
        elseif j == 7
            ext_chl = ext_chl - NN
            if ext_chl < 0.0
              ext_chl = 0.0
            end
            int_chl = int_chl + 1
            asum = asum - a[7] - a[8]
            a[7] = D_1 * ext_chl * NN /KK
            a[8] = K_5 * ((np / V) / (1 + K_6 / (int_chl / V)))
            asum = asum + a[7] + a[8]
        elseif j == 8
            int_chl = int_chl - 1
            asum = asum - a[8]
            a[8] = K_5 * ((np / V) / (1 + K_6 / (int_chl / V)))
            asum = asum + a[8]
        elseif j == 9
            np = np - 1
            asum = asum - a[8] - a[9]
            a[8] = K_5 * ((np / V) / (1 + K_6 / (int_chl / V)))
            a[9] = d2 * np
            asum = asum + a[8] + a[9]
        elseif j == 10
            rp = rp - 1
            asum = asum - a[10]
            a[10] = d2 * rp
            asum = asum + a[10]
        end
    end
end

# define code which ties SSA to model of fixed number of cells growing and dividing using logistic growth
function constitutive_code(OD::Float64, KK::Int64,init_protein::Int64,init_mRNA::Int64,L::Float64,maxcell::Int64,timestep::Float64,maxtime::Float64,Ni::Float64,Nd::Float64,Nf::Float64,b1::Float64,b2::Float64,q1::Float64,q2::Float64,d1::Float64,d2::Float64,mu0::Float64,D_1::Float64,K_M::Float64,K_5::Float64,K_6::Float64, chl_pulse::Array{Float64,1},max_chl_cell::Float64,delay::Float64,dd::Float64)

  output_V = zeros(maxcell,Int(maxtime/timestep))
  output_N = zeros(maxcell,Int(maxtime/timestep))
  output_R = zeros(maxcell,Int(maxtime/timestep))
  output_NP = zeros(maxcell,Int(maxtime/timestep))
  output_RP = zeros(maxcell,Int(maxtime/timestep))
  output_mu = zeros(maxcell,Int(maxtime/timestep))
  output_ext_chl = zeros(maxcell,Int(maxtime/timestep))
  output_int_chl = zeros(maxcell,Int(maxtime/timestep))
  init_cell_count = maxcell
  final_length = zeros(maxcell*2)
  n= zeros(Int64,maxcell*2)
  r = zeros(Int64,maxcell*2)
  np= zeros(Int64,maxcell*2)
  rp = zeros(Int64,maxcell*2)
  delays = zeros(maxcell)
  V = zeros(maxcell*2)
  mu = zeros(maxcell*2)
  division_dist = Normal(0.5,Nd)
  delay_dist = Normal(delay,dd)
  ext_chl = zeros(Int, maxcell)
  int_chl = zeros(Int, maxcell * 2)
  number_divisions = round(OD * KK)
  NN = ones(Int(maxtime / timestep)) .* number_divisions
  count = 1
  for i = 1:init_cell_count
    V[i] = L+ Ni*L*rand()
    final_length[i] = 1+ L  + Nf*randn()
    n[i] = Int(round(init_mRNA*rand()))
    r[i] = n[i]
    np[i] = Int(round(init_protein*V[i]))
    rp[i] = Int(round(init_protein*V[i]))
    mu[i] = 0.0
    delays[i] = rand(delay_dist)
  end
  count_error = 0
  extra_cells = 0
  pop2Sample = collect(1:maxcell)
  output_V[:,count] = V[1:maxcell]
  output_N[:,count] = n[1:maxcell]
  output_R[:,count] = r[1:maxcell]
  output_NP[:,count] = np[1:maxcell]./V[1:maxcell]
  output_RP[:,count] = rp[1:maxcell]./V[1:maxcell]
  output_mu[:,count] = mu[1:maxcell]
  output_int_chl[:,count] = int_chl[1:maxcell]
  s = chl_pulse
  for ts = timestep:timestep:maxtime-timestep
    i1 = 0
    while i1 < maxcell
      i1 = i1 + 1
      ext_chl[i1] = ext_chl[i1] + Int(round(s[count] .* max_chl_cell))

      nn, nnp, rr, rrp, intchl,extchl = ssa_antibiotic(b1, q1, mu[i1], b2, q2, D_1, K_5, K_6, V[i1], timestep, d1, d2, n[i1], np[i1], r[i1], rp[i1], Int(NN[count]), ext_chl[i1], int_chl[i1],KK)      
      n[i1] = nn
      np[i1] = nnp
      r[i1] = rr
      rp[i1] = rrp
      int_chl[i1] = intchl
      ext_chl[i1] = extchl


      mu[i1] = mu0 * (1 / (1 + (int_chl[i1]) / (V[i1] * K_M))) * (NN[count]/KK) * (1.0 - (NN[count] / KK))

      V[i1] = V[i1]*exp(mu[i1]*timestep)

      if V[i1] > final_length[i1] && number_divisions < KK
        
        number_divisions = number_divisions + 1

        extra_cells = extra_cells + 1
        new_cell = extra_cells + maxcell
        R = rand(division_dist)
        V[new_cell] = V[i1]*R
        V[i1] = V[i1]*(1.0-R)

        temp_n = n[i1]
        temp_r = r[i1]
        daughter_n = rand(Binomial(Int(temp_n),R))
        daughter_r = rand(Binomial(Int(temp_r),R))
        n[i1] = daughter_n
        r[i1] = daughter_r
        n[new_cell] = temp_n - daughter_n
        r[new_cell] = temp_r - daughter_r

        temp_np = np[i1]
        temp_rp = rp[i1]
        daughter_np = rand(Binomial(Int(temp_np),R))
        daughter_rp = rand(Binomial(Int(temp_rp),R))
        np[i1] = daughter_np
        rp[i1] = daughter_rp
        np[new_cell] = temp_np - daughter_np
        rp[new_cell] = temp_rp - daughter_rp

        temp_int_chl = int_chl[i1]
        daughter_int_chl = rand(Binomial(Int(temp_int_chl),R))
        int_chl[i1] = daughter_int_chl
        int_chl[new_cell] = temp_int_chl - daughter_int_chl

        mu[i1] = mu0 * (1 / (1 + (int_chl[i1]) / (V[i1] * K_M))) * (NN[count] / KK) * (1.0 - (NN[count] / KK))
        mu[new_cell] = mu0 * (1 / (1 + (int_chl[new_cell]) / (V[new_cell] * K_M))) * (NN[count] / KK) * (1.0 - (NN[count] / KK))

        final_length[i1] = 1+V[i1] + Nf*randn()
        final_length[new_cell] = 1+V[new_cell] + Nf*randn()
      end
    end
    if extra_cells > 0
      remov_index = shuffle(pop2Sample)
      remov = remov_index[1:extra_cells]
      for k = 1:extra_cells
        V[remov[k]] = V[maxcell+ k]
        n[remov[k]] = n[maxcell+ k]
        r[remov[k]] = r[maxcell+ k]
        np[remov[k]] = np[maxcell+ k]
        rp[remov[k]] = rp[maxcell+ k]
        mu[remov[k]] = mu[maxcell+ k]
        int_chl[remov[k]] = int_chl[maxcell+ k]
        final_length[remov[k]] = final_length[maxcell+ k]
      end
      extra_cells = 0
    end
    count = count + 1
    output_V[:,count] = V[1:maxcell]
    output_N[:,count] = n[1:maxcell]
    output_R[:,count] = r[1:maxcell]
    output_NP[:,count] = np[1:maxcell]./V[1:maxcell]
    output_RP[:,count] = rp[1:maxcell]./V[1:maxcell]
    output_int_chl[:,count] = int_chl[1:maxcell]
    output_ext_chl[:,count] = ext_chl[1:maxcell]
    output_mu[:,count] = mu[1:maxcell]
    NN[count] = number_divisions
  end
  return output_V,output_N,output_R,output_NP,output_RP,count_error,output_mu,number_divisions,output_int_chl, NN,output_ext_chl
end

