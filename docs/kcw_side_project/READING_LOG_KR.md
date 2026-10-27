# K_cw 사이드 프로젝트 — 읽은 자료와 찾은 기록

받은 논문·SI·기록과 거기서 얻은 것을 모두 적습니다. 원문 PDF와 DOCX는 저작권 때문에 커밋하지 않습니다
(`docs/literature_db/raw_si/README.md`). 대신 추출값을 `raw_si/`에 두고, 원본은 아래 SHA-256으로 확인합니다.
원문 보관처: Google Drive 폴더 `pfas-rice-model`(사용자 소유, 2026-10-27 생성). 파일 목록과 SHA-256은 그 폴더의
`kcw_papers_file_index.csv`에 있습니다. 커넥터는 대용량 파일을 올릴 수 없어서 원문 파일은 사용자가 직접 올립니다.

## 1. 원문과 SI (2026-10-27)

| 자료 | 받은 경로 | 읽은 범위 | 얻은 것 → 저장 위치 |
|---|---|---|---|
| Mel 2024, J. Hazard. Mater. `10.1016/j.jhazmat.2024.136016` | 사용자 업로드(papers_5.zip) + SI(papers_SI.zip) | 본문 전체, SI(Table S1은 그림이라 값을 읽지 못함, Fig. S1 Langmuir R²) | 리그닌 log Kd 4종, Freundlich n·KF, 흡착·탈착 % → `raw_si/mel2024_lignin_kd.csv` |
| Guo 2025, Nat. Commun. `10.1038/s41467-025-65191-3` | 오픈 액세스(nature.com) | 본문의 관련 절, Source Data(MOESM7) Fig. 3b | DFT는 PFOA 하나뿐. 세포 이하 분획 10종 → `raw_si/guo2025_fig3b_subcellular.csv`. SI PDF(MOESM1) Table 11은 아직 안 읽음 |
| Muensterman 2026, ES&T `10.1021/acs.est.5c11716` + 정정 `10.1021/acs.est.6c07668` | 업로드 | 본문, 정정문, SI Table S10–S13 | S12 조각 기여 → `raw_si/muensterman2026_tableS12_units.csv`. S13의 PFCA·PFSA·에테르 39종 → `raw_si/muensterman2026_tableS13_kiam.csv` |
| Wang 2020, Sci. Total Environ. `10.1016/j.scitotenv.2020.139383` | 업로드 | 본문, SI 표 전부 | 분획 농도·비율·소실, 문헌표 → `raw_si/wang2020_tables10_conc.csv`, `_tables11_proportion.csv`, `_tables13_dissipation.csv`, `_tables15_literature.csv` |
| Liu 2023, ES&T `10.1021/acs.est.3c00504` | 업로드 | 본문, SI | 뿌리 PFAS의 81.5–91.2 %가 표면 흡착 분획(`K_surf` 쪽). 분획 비율은 본문 Fig. 3a에만 있음(디지털화 필요). SI에는 분획 표가 없음 |
| Yu 2021, ES&T `10.1021/acs.est.1c01175` | 업로드 | 본문, SI | 상추 PFOS 세포벽(F1) 46–62 %. 값은 Fig. 1과 SI Fig. S6(디지털화 필요) |
| Mei 2021, Environ. Int. `10.1016/j.envint.2021.106642` | 업로드(ScienceDirect는 막힘) | 본문 | 리뷰. 세포벽 성분(셀룰로오스·펙틴·리그닌)의 Kd는 없음. K_cw_poly의 출처가 아님 |

SI 원본(사용자 manifest `si_manifest.csv`):

| 파일 | 원래 이름 | SHA-256 |
|---|---|---|
| Mel2024_SI1.docx | 1-s2.0-S0304389424025950-mmc1.docx | e42484ca6467230c138c3ab27b1850cf02ffd5266b7e28b18e037276f2c87701 |
| Muensterman2025_SI1.pdf | es5c11716_si_001.pdf (figshare 60722189) | ca08dfc15ceeddebe249d3bdb259403c8347d8460d418afc75e5e7802ce1fa0d |
| Liu2023_SI1.pdf | es3c00504_si_001.pdf (figshare 41001883) | 1408f968ff21e0ae18159447354aa0b1ece6d78dd913364498f9d68e27cf004b |
| Wang2020_SI1.docx | 1-s2.0-S0048969720329004-mmc1.docx | 4cc0165b75a4ac0a1b747d04d434828d8cbc56dcca829128ede27a97082cb985 |
| Yu2021_SI1.pdf | es1c01175_si_001.pdf (figshare 28580577) | 157f0b57ec1bc07fce21c502030787f31c6974f959bcc87031847c823180b4af |

## 2. 읽고 바뀐 판단
- **Mel 2024가 GAP A 리그닌 4점과 맞습니다.** CSV의 low–high는 log Kd ± SE입니다.
- **PFSA 기울기가 다릅니다.** Mel 두 점에서 PFCA는 0.085 log/CF₂, PFSA는 0.028 log/CF₂입니다. 사다리는 둘 다 0.10을 씁니다.
- **설폰산 오프셋은 사슬길이에 따라 줄어듭니다.** CF₂를 맞추면 C4 +0.47, C8 +0.24 log입니다.
- **에테르 결합**(Muensterman Table S12): 지질 친화성(CHI log D)은 +0.10…+0.49로 오르지만, 막 분배(log K_IAM)는 −0.06…−0.31로 내려갑니다.
  막 분배의 부호는 Chen 2025의 K_MW 오프셋(−0.49/O)과 같습니다. 크기는 더 작습니다.
  (2026-10-27 오전의 계획서에 "부호가 반대"라고 쓴 것은 본문만 보고 CHI log D를 읽은 것이었고, SI를 보고 고쳤습니다.)
- **Guo 2025 Supplementary Table 11**(SI p. S42)은 *O. rosea* 뿌리·잎의 펙틴과 헤미셀룰로오스 단당 조성만 담고 있습니다. GAP A 문서
  "Experimental design"과 H8 handoff가 말한 "Guo Table 11 벼판(rice panel)"은 이 논문에 없습니다.
- **디지털화**(2026-10-27): Liu 2023 Fig. 3a(밀 뿌리, 4종 ± 부식산)와 Yu 2021 Fig. 1b(상추 뿌리 PFOS)를 300 dpi 렌더링의 색 구간으로 읽었습니다.
  본문 범위와 맞습니다(Liu F1 47.8–59.0 % 대 본문 47.74–58.96 %, Yu F1 46.4–62.0 % 대 본문 46–62 %).
  Liu 뿌리 세포벽 비율: PFHxS 0.50, PFOS 0.56, PFOA 0.57, 6:2 Cl-PFESA 0.59.
- **Muensterman Table S13의 HFPO-DA log K_IAM은 2.0**으로 PFHxA(2.0)와 같고, ADONA는 2.3입니다. 막 분배 쪽 비교값입니다(세포벽 값은 아님).

## 3. 찾지 못한 기록: "Prompt2 handoff"
- 찾는 것: `params/Kcw_Klignin_params_v2.csv`와 함께 왔다는 handoff. K_cw_poly 절대 크기(PFOA 4.0), 설폰산 오프셋 0.40, 기울기 0.10의 출처.
- 찾은 곳:
  - 저장소와 git 이력(`git log --all -S`)
  - Google Drive, Notion
  - 패키지를 받은 Claude Code 세션 `session_01CPrFZKenAkXeWM12MM4PZX`의 기록. 2026-06-10 04:49부터 06-11 03:06까지 사용자·도우미 메시지를 전부 넘겨 보았습니다.
    16:53:36 업로드에는 `pfas_rice_model.zip`, `README 1.md`, `parameters.json`만 있었고, 압축 파일에도 handoff는 없었습니다.
  - 사용자 확인(2026-10-27): 구하지 못함.
- 판단: `docs/H8_handoff_S6_final.md` 136행의 "본 세션"은 H8을 쓴 바깥 대화(S6)입니다. handoff는 거기에만 있었던 것으로 보입니다.
- 단서였던 문헌: Mei 2021(성분 Kd 없음), Wang 2020(성분 Kd 없음), "cellulose/lignin biosorbent(선택, 미정)"(구체적 문헌 없음).
- 결론: K_cw_poly는 출처를 추적할 수 없는 값입니다(`PLAN_KR.md` §2.4).
