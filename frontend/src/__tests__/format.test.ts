import { describe, expect, it } from "vitest";
import fixtures from "./format_fixtures.json";
import { fmtG, pyFormat } from "../format";

describe("fmtG matches Python format(v, '.pg')", () => {
  it("all fixture values", () => {
    for (const [v, p, want] of fixtures as [number, number, string][]) {
      expect(`${v}:.${p}g -> ${fmtG(v, p)}`).toBe(`${v}:.${p}g -> ${want}`);
    }
  });
});

describe("pyFormat", () => {
  it("formats named fields with specs", () => {
    expect(pyFormat("95% range {lo:.2g}–{hi:.2g}", { lo: 0.01234, hi: 123.4 })).toBe("95% range 0.012–1.2e+02");
    expect(pyFormat("{name}: you {meas:.3g}", { name: "root", meas: 0.47898 })).toBe("root: you 0.479");
    expect(pyFormat("{missing}", {})).toBe("{missing}");
  });
});

import ffix from "./fmtf_fixtures.json";
import { fmtF } from "../format";

describe("fmtF matches Python format(v, '.df')", () => {
  it("all fixture values (incl. exact .5 ties -> half-even)", () => {
    for (const [v, d, want] of ffix as [number, number, string][]) {
      expect(`${v}:.${d}f -> ${fmtF(v, d)}`).toBe(`${v}:.${d}f -> ${want}`);
    }
  });
});
