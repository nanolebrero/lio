{
  scalar_type F_mU[2];
  {
    scalar_type U =
        (PmC[0] * PmC[0] + PmC[1] * PmC[1] + PmC[2] * PmC[2]) * (ai + aj);
    // TODO (maybe): test out storing F(m,U) values in texture and doing a
    // texture fetch here rather than the function calculation
    lio_gamma<scalar_type, 1>(F_mU, U);
  }
  {
    // START INDEX i1=0, CENTER 1
    {
      scalar_type p1s_0 = PmA[0] * F_mU[0] - PmC[0] * F_mU[1];
      my_fock[0] += (double)(clatom_charge_sh[j] * p1s_0);
    }
    // START INDEX i1=1, CENTER 1
    {
      scalar_type p1s_0 = PmA[1] * F_mU[0] - PmC[1] * F_mU[1];
      my_fock[1] += (double)(clatom_charge_sh[j] * p1s_0);
    }
    // START INDEX i1=2, CENTER 1
    {
      scalar_type p1s_0 = PmA[2] * F_mU[0] - PmC[2] * F_mU[1];
      my_fock[2] += (double)(clatom_charge_sh[j] * p1s_0);
    }
  }
}
