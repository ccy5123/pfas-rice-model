"""
Plant-air exchange: volatilisation and gaseous uptake
=====================================================

The core ODE (`pfas_rice_plant_module_4pool_surf`) has no air terms at all.  For
PFAS that is a deliberate and correct simplification -- a perfluorinated anion has
`K_AW ~ 0`, so there is nothing to volatilise (CLAUDE.md section 2 turns air
exchange off explicitly).  But the same ODE now also carries the NEUTRAL organic
path (`neutral_dpu`), which is meant to cover ordinary organics, and many of those
ARE volatile.  Without air terms every such compound is an upper bound by
construction: `neutral_dpu_validation` section 3 shows the leaf is an unbounded
terminal accumulator whose only sinks are metabolism and volatilisation, and until
now only one of the two existed.

This module supplies the missing half.  It implements the equation set already
written out in `docs/dpu_model_summary_corrected.tex` section `sec:permeability`:

    eq:Pc     cuticle permeability          -> p_cuticle
    eq:Pair   air boundary-layer conductance-> p_air_boundary
    eq:Paqua  aqueous/apoplast permeability -> p_aqueous
    eq:Pctot  the three in SERIES           -> p_cuticle_total
    eq:Ps     stomatal conductance          -> p_stomata
    eq:Csat   saturation water vapour       -> c_h2o_sat
    eq:Pp     cuticle + stomata in PARALLEL -> AirExchange.p_plant
    eq:Kpa    plant-air partition           -> AirExchange.k_pa
    eq:Qvol   volatilisation flux           -> AirExchange.flux (loss term)
    eq:Qgas   gaseous uptake                -> AirExchange.flux (gain term)

Three things the derivation flags, and how they are handled here
----------------------------------------------------------------
1. **Units.** The numerical constants in the permeability correlations assume SI
   BASE units (m, s, g/mol), while this repo works in day / L / kg / ug.  Every
   public function here returns **m/day**, converted once at the source, and the
   flux is returned in the ODE's own **ug/(kg day)**.  The m^3 -> L conversion in
   the flux is the `M3_TO_L` factor; getting it wrong is a 1000x error, so it is
   isolated and tested.
2. **The `1/(1-phi)` singularity** in eq:Ps blows up as relative humidity -> 1.
   The tex notes the product stays finite because the transpiration stream -> 0 in
   that limit, but the factor itself must be guarded: `rh` is capped at `RH_MAX`.
3. **It must vanish for PFAS.**  This is structural rather than a special case:
   `P_air` is proportional to `K_AW` and enters eq:Pctot as a SERIES resistance, so
   `K_AW = 0` makes the whole cuticular path zero; `P_S` is likewise proportional
   to `K_AW`.  Both pathways therefore vanish, and the volatilisation term carries
   a further factor `K_AW` through `1/K_PA`.  On top of that the core keeps
   `air=None` as its default, so the terms are not merely zero but never evaluated
   -- every existing PFAS number, and `reproduce_demo`'s RMSE 0.029, is
   bit-identical.  (Same discipline as `k_seq=0` in the two-pool merge.)

Modelling assumptions, quoted from the derivation
-------------------------------------------------
* **No volatilisation from the roots** (below-ground tissue).
* **Stem: cuticle only** (no stomata).
* **Leaf and fruit/grain: cuticle and stomata in parallel.**

Scope and honesty
-----------------
* **Specific surface area `S` is the weak input, not the equations.**  The flux
  scales linearly with the tissue's specific surface `S` [m^2/kg], and `S` entered
  this repo as a leaf/grain RATIO for splitting the xylem stream -- only the ratio
  was ever load-bearing, so the absolute values (leaf 20, grain 2 m^2/kg) have
  never been calibrated as areas.  Using them as absolute areas is a new demand on
  them.  `AirExchange.S` therefore exists to override them with measured values,
  and the shipped rice compartments give the **stem S = 0**, which makes the stem
  term inert until a real stem surface area is supplied.  Treat the ABSOLUTE
  magnitude of a volatilisation loss from this module as order-of-magnitude until
  the areas are measured; the DIRECTION and the Kow/K_AW dependence are the
  derivation's.
* **Particle deposition (eq:Qdep) is NOT implemented.**  It is a separate
  atmospheric-deposition pathway needing a particle-bound air concentration and a
  deposition velocity, not part of the plant-air equilibrium exchange, and the
  `f_particle` fraction here only excludes the particle-bound share from the
  gaseous uptake, as eq:Qgas specifies.
* Only the 4-compartment core (`RiceUptakeModel`) carries the hook.  The
  multi-segment shoot (`NStemLeafModel`) does not -- that is a PFAS-side model and
  PFAS have no air exchange.

References: Trapp & Matthies (1995) ES&T 29:2333-2338 for the cuticle/stomata
formulation; the equation labels above point at this repo's own derivation, which
is where the constants used here are stated.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from pfas_rice_plant_module_4pool_surf import ROOT, STEM, LEAF, FRUIT

# --- unit conversions (the SI -> repo bridge; see note 1 in the header) --------
SEC_PER_DAY = 86400.0     # s/day: the correlations are per SECOND
M3_TO_L = 1000.0          # L/m^3: areas are m^2 and permeabilities m/day, but
                          # concentrations are per LITRE

# --- physical constants -------------------------------------------------------
RHO_WATER = 1000.0        # kg/m^3
R_WATER_VAPOUR = 461.9    # J/(kg K), specific gas constant of water vapour (eq:Csat)
D_O2_WATER = 2.1e-9       # m^2/s, diffusivity of O2 in water (reference solute of eq:Paqua)
MW_O2 = 32.0              # g/mol
MW_H2O = 18.0             # g/mol
RH_MAX = 0.999            # cap on relative humidity (guards the 1/(1-phi) pole)

# --- eq:Pc coefficients -------------------------------------------------------
CUTICLE_LOGKOW_SLOPE = 0.704
CUTICLE_INTERCEPT = -11.2


def p_cuticle(log_kow: float) -> float:
    """Cuticle permeability (eq:Pc): `10^(0.704 log Kow - 11.2)` [m/day].

    The correlation is in m/s; converted here once.  Note it does NOT involve
    K_AW -- the cuticle is a lipid phase, so a non-volatile compound still has a
    finite cuticular permeability.  What makes the whole cuticular PATH vanish for
    a non-volatile compound is the air boundary layer in series with it (eq:Pair).
    """
    return 10.0 ** (CUTICLE_LOGKOW_SLOPE * float(log_kow) + CUTICLE_INTERCEPT) * SEC_PER_DAY


def p_air_boundary(K_AW: float, MW: float) -> float:
    """Air boundary-layer conductance (eq:Pair): `sqrt(300) K_AW / (200 sqrt(MW))` [m/day].

    Proportional to K_AW, which is what drives the whole air pathway to zero for a
    PFAS-like compound: it sits in SERIES with the cuticle (eq:Pctot), so a zero
    here zeroes the cuticular path regardless of how permeable the cuticle is.
    """
    if MW <= 0:
        raise ValueError(f"molar mass must be positive, got {MW}")
    return np.sqrt(300.0) * float(K_AW) / (200.0 * np.sqrt(float(MW))) * SEC_PER_DAY


def p_aqueous(MW: float, z_path: float = 1e-3, D_O2: float = D_O2_WATER) -> float:
    """Aqueous (apoplast) permeability (eq:Paqua): `(D_O2/z) (MW/32)^(-1/2)` [m/day].

    Oxygen is the reference solute; the molar-mass ratio scales its diffusivity to
    the compound of interest.  `z_path` is the diffusion path length [m] -- the
    thickness of the stagnant aqueous layer -- and is the one free geometric input
    (default 1 mm).  It is rarely limiting: for a lipophilic compound the cuticle
    is orders of magnitude slower.
    """
    if z_path <= 0:
        raise ValueError(f"diffusion path length must be positive, got {z_path}")
    return (D_O2 / float(z_path)) * (float(MW) / MW_O2) ** -0.5 * SEC_PER_DAY


def _series(*perms: float) -> float:
    """Combine permeabilities acting in SERIES: `(sum 1/P)^-1`, zero-safe.

    A zero member is an infinite resistance, so the series total is zero -- which
    is exactly the K_AW = 0 case (eq:Pair -> 0).  Computing it as `1/sum(1/P)`
    would divide by zero, so it is short-circuited here rather than relying on
    IEEE inf arithmetic.
    """
    ps = [float(p) for p in perms]
    if any(p <= 0.0 for p in ps):
        return 0.0
    return 1.0 / sum(1.0 / p for p in ps)


def p_cuticle_total(log_kow: float, K_AW: float, MW: float,
                    z_path: float = 1e-3, D_O2: float = D_O2_WATER) -> float:
    """Total cuticular permeability (eq:Pctot): cuticle, air layer and aqueous
    layer in SERIES [m/day].  Zero when `K_AW = 0`."""
    return _series(p_cuticle(log_kow),
                   p_air_boundary(K_AW, MW),
                   p_aqueous(MW, z_path, D_O2))


def c_h2o_sat(T: float = 298.15) -> float:
    """Saturation water-vapour concentration (eq:Csat) [kg/m^3].

    `p_H2O = 610.7 * 10^(7.5 T_C / (237 + T_C))` [Pa] via the ideal gas law with
    the specific gas constant of water vapour.  ~0.023 kg/m^3 at 25 C.
    """
    T = float(T)
    T_C = T - 273.15
    p = 610.7 * 10.0 ** (7.5 * T_C / (237.0 + T_C))
    return p / (R_WATER_VAPOUR * T)


def p_stomata(K_AW: float, MW: float, q_per_area: float,
              rh: float = 0.7, T: float = 298.15) -> float:
    """Stomatal conductance (eq:Ps) [m/day].

        P_S = rho_water / ((1 - phi) C_H2O,sat) * (MW/18)^(-1/2) * K_AW * (Q_XYL/A)

    The stomatal pathway is *tied to the transpiration stream*: a compound leaves
    through the stomata in proportion to how much water is leaving, scaled by its
    volatility relative to water.  `q_per_area` is the organ's share of the xylem
    flow per unit surface area [m/day] -- i.e. the transpiration velocity.

    `rh` is capped at RH_MAX: the `1/(1-phi)` factor is singular at saturation.
    The tex notes the product stays finite there because `Q_TP -> 0` as well, so
    the cap only guards the arithmetic, it does not change the physics.
    """
    phi = min(float(rh), RH_MAX)
    if phi < 0.0:
        raise ValueError(f"relative humidity must be in [0,1), got {rh}")
    return (RHO_WATER / ((1.0 - phi) * c_h2o_sat(T))
            * (float(MW) / MW_H2O) ** -0.5 * float(K_AW) * float(q_per_area))


@dataclass
class AirExchange:
    """Plant-air exchange for one compound: volatilisation + gaseous uptake.

    Attach to `RiceUptakeModel(air=...)`; leaving it `None` (the default) skips the
    terms entirely, so every PFAS result is untouched.

    Parameters
    ----------
    K_AW : air-water partition coefficient [-] (dimensionless Henry's law
        constant).  `0` makes every term identically zero, by construction.
    MW : molar mass [g/mol].  Required -- eq:Pair, eq:Paqua and eq:Ps all need it.
    log_kow : octanol-water partition coefficient, log10 [-]  (eq:Pc).
    C_air : ambient GAS-phase concentration [ug/m^3].  Default 0 = clean air, i.e.
        volatilisation only.  Note the unit: air concentrations are conventionally
        per cubic metre while the rest of the model is per litre.
    f_particle : fraction of the air concentration bound to particles [-]
        (eq:Qgas).  Only the gaseous `(1 - f_particle)` share is taken up here;
        particle deposition (eq:Qdep) is not implemented -- see the header.
    rh : relative humidity [-], capped at RH_MAX.
    T : temperature [K] (eq:Csat).
    z_path : aqueous boundary-layer thickness [m] (eq:Paqua).
    D_O2 : diffusivity of oxygen in water [m^2/s] (eq:Paqua).
    S : optional per-organ specific surface area [m^2/kg fw], overriding the
        compartments' own `S`.  Keys: 'stem', 'leaf', 'grain'.  READ THE HEADER --
        the compartments' `S` was introduced as a leaf/grain ratio, not as a
        calibrated absolute area, and the shipped stem `S` is 0.
    """
    K_AW: float
    MW: float
    log_kow: float
    C_air: float = 0.0
    f_particle: float = 0.0
    rh: float = 0.7
    T: float = 298.15
    z_path: float = 1e-3
    D_O2: float = D_O2_WATER
    S: dict | None = None
    # organs that exchange with air, and whether they have stomata (eq:Pp and the
    # derivation's "important modelling assumptions"): the root never does.
    _STOMATAL = (LEAF, FRUIT)
    _EXCHANGING = (STEM, LEAF, FRUIT)

    def __post_init__(self):
        if not np.isfinite(self.MW) or self.MW <= 0:
            raise ValueError(
                "AirExchange needs a molar mass [g/mol]: eq:Pair, eq:Paqua and "
                f"eq:Ps all scale with it (got MW={self.MW}). Supply "
                "NeutralCompound(MW=...) before enabling air exchange.")
        if self.K_AW < 0:
            raise ValueError(f"K_AW must be >= 0, got {self.K_AW}")

    # -- permeabilities --------------------------------------------------------
    def p_cuticular(self) -> float:
        """eq:Pctot for this compound [m/day] (organ-independent)."""
        return p_cuticle_total(self.log_kow, self.K_AW, self.MW,
                               self.z_path, self.D_O2)

    def p_plant(self, q_per_area: float = 0.0, stomata: bool = False) -> float:
        """Total plant-air permeability (eq:Pp) [m/day].

        Cuticle and stomata act in PARALLEL, so they ADD (unlike the three
        resistances inside eq:Pctot, which are in series).
        """
        p = self.p_cuticular()
        if stomata:
            p += p_stomata(self.K_AW, self.MW, q_per_area, self.rh, self.T)
        return p

    def k_pa(self, K_PW: float, rho: float = 1.0) -> float:
        """Plant-air partition coefficient (eq:Kpa): `K_PW * rho / K_AW` [-].

        Reporting only.  The flux does NOT go through this function, because
        `rho` cancels analytically against eq:Qvol's `A*rho/M` prefactor (see
        `flux`) and because K_PA is infinite at K_AW = 0, which is precisely the
        case the flux has to evaluate without blowing up.
        """
        if self.K_AW <= 0:
            return float("inf")
        return float(K_PW) * float(rho) / float(self.K_AW)

    # -- fluxes ----------------------------------------------------------------
    def specific_area(self, comps, idx: int) -> float:
        """Specific surface area [m^2/kg fw] of compartment `idx`."""
        if self.S is not None:
            return float(self.S.get(comps[idx].name, 0.0))
        return float(comps[idx].S)

    def flux(self, C, M, K_PW, comps, Qtp: float = 0.0,
             xyl_share=(0.0, 0.0)) -> np.ndarray:
        """Net air exchange per compartment [ug/(kg day)]: gaseous uptake MINUS
        volatilisation, ready to add to the ODE's `dC/dt`.

        Parameters
        ----------
        C : tissue concentrations [ug/kg], length 4.
        M : tissue fresh masses [kg], length 4.
        K_PW : plant-water partition per tissue [L/kg] -- the `binding_factors`
            of the run.  For a neutral compound that IS the Briggs K_PW; the
            equations are indifferent to which model produced the number.
        Qtp : transpiration stream [L/day].
        xyl_share : (leaf, grain) fractions of `Qtp` reaching each organ, i.e. the
            core's surface-area split (and, for the grain, its formation gate).

        Units, derived once so the 1000 is not a magic number:
            volatilisation   (A rho / M) * (P C / K_PA)                  (eq:Qvol)
                           = S rho P C K_AW / (K_PW rho)                 (eq:Kpa)
                           = S P C K_AW / K_PW           <- rho cancels exactly
              [m^2/kg][m/day][ug/kg][-]/[L/kg] = m^3/(kg day) * ug/L
              -> multiply by M3_TO_L to reach ug/(kg day).
            gaseous uptake   (1 - f_p) C_A (A/M) P                       (eq:Qgas)
              [ug/m^3][m^2/kg][m/day] = ug/(kg day) directly -- no factor, because
              C_air is per CUBIC METRE.

        The two balance (net flux zero) when the air is in equilibrium with the
        tissue, `C_air = M3_TO_L * K_AW * C / K_PW`, which is Henry's law on the
        tissue's free aqueous concentration `C/K_PW`.  `test_plant_air` pins that.
        """
        out = np.zeros(4)
        if self.K_AW <= 0.0:
            return out                      # no air pathway at all (PFAS)
        C = np.asarray(C, dtype=float)
        M = np.asarray(M, dtype=float)
        shares = {LEAF: float(xyl_share[0]), FRUIT: float(xyl_share[1])}
        for idx in self._EXCHANGING:
            S_k = self.specific_area(comps, idx)
            if S_k <= 0.0:
                continue                    # no area -> no exchange (stem default)
            stomata = idx in self._STOMATAL
            q_per_area = 0.0
            if stomata:
                A = S_k * max(float(M[idx]), 0.0)      # [m^2]
                if A > 0.0:
                    # organ's share of the transpiration stream per unit area:
                    # L/day -> m^3/day -> divided by m^2 gives m/day
                    q_per_area = (shares[idx] * float(Qtp) / M3_TO_L) / A
            P = self.p_plant(q_per_area, stomata=stomata)
            if P <= 0.0:
                continue
            gain = (1.0 - self.f_particle) * float(self.C_air) * S_k * P
            loss = M3_TO_L * S_k * P * C[idx] * self.K_AW / float(K_PW[idx])
            out[idx] = gain - loss
        return out

    # -- reporting -------------------------------------------------------------
    def summary(self, comps, K_PW, M=None, Qtp: float = 0.0,
                xyl_share=(1.0, 0.0)) -> dict:
        """Per-organ permeabilities and the implied volatilisation half-life.

        The half-life is the diagnostic that makes the term interpretable: it is
        the time for the organ to lose half its burden to clean air at fixed mass,
        `ln2 / (M3_TO_L * S * P * K_AW / K_PW)`, i.e. directly comparable to the
        metabolic half-life the neutral path already scans.

        It is a SNAPSHOT at the `(M, Qtp)` passed in, not a season average: the
        stomatal pathway follows the transpiration stream, so the same compound
        volatilises faster at peak transpiration than at maturity. Read it as the
        order of magnitude of the sink, not as an integrated loss.
        """
        M = np.ones(4) if M is None else np.asarray(M, dtype=float)
        shares = {LEAF: float(xyl_share[0]), FRUIT: float(xyl_share[1])}
        rows = {}
        for idx in self._EXCHANGING:
            S_k = self.specific_area(comps, idx)
            stomata = idx in self._STOMATAL
            q_per_area = 0.0
            if stomata and S_k > 0 and M[idx] > 0:
                q_per_area = (shares[idx] * float(Qtp) / M3_TO_L) / (S_k * M[idx])
            P = self.p_plant(q_per_area, stomata=stomata)
            k = M3_TO_L * S_k * P * self.K_AW / float(K_PW[idx])
            rows[comps[idx].name] = dict(
                S=S_k, P_cuticular=self.p_cuticular(),
                P_stomatal=(p_stomata(self.K_AW, self.MW, q_per_area, self.rh, self.T)
                            if stomata else 0.0),
                P_plant=P, K_PA=self.k_pa(K_PW[idx]),
                rate=k, half_life=(np.log(2.0) / k if k > 0 else float("inf")))
        return rows


def _demo():
    """Permeabilities and volatilisation half-lives across volatility."""
    import neutral_dpu as ND

    comps = ND.rice_compartments()
    print(f"{'compound':>16}{'logKow':>8}{'K_AW':>10}{'P_cut':>11}{'P_stom':>11}"
          f"{'leaf t1/2':>11}")
    print(f"{'':>16}{'':>8}{'':>10}{'[m/d]':>11}{'[m/d]':>11}{'[d]':>11}")
    # a volatility ladder at fixed lipophilicity, plus PFOA for the null case
    for name, lk, kaw, mw in (("PFOA (anion)", 4.8, 0.0, 414.0),
                              ("carbamazepine", 2.45, 4.6e-9, 236.3),
                              ("chlorpyrifos", 4.01, 4.8e-4, 350.6),
                              ("trichloroethene", 2.42, 0.40, 131.4)):
        air = AirExchange(K_AW=kaw, MW=mw, log_kow=lk)
        c = ND.neutral_compound(ND.NeutralCompound(name, lk))
        from pfas_rice_plant_module_4pool_surf import binding_factors
        K_PW = binding_factors(comps, c)
        s = air.summary(comps, K_PW, M=np.array([0.03, 0.04, 0.05, 0.025]),
                        Qtp=0.1, xyl_share=(0.9, 0.1))["leaf"]
        hl = f"{s['half_life']:.3g}" if np.isfinite(s["half_life"]) else "inf"
        print(f"{name:>16}{lk:>8.2f}{kaw:>10.1e}{s['P_cuticular']:>11.2e}"
              f"{s['P_stomatal']:>11.2e}{hl:>11}")
    print("\nK_AW = 0 (PFAS) gives an infinite half-life: the air pathway is")
    print("structurally absent, not numerically small.")


if __name__ == "__main__":
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    _demo()
