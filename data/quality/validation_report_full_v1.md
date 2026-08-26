# 8-oxoG 염기쌍 내부 변수 좌표 재계산 검증 보고서 v1

> **후속 상태 메모(2026-08-26):** 이 문서는 2026-08-04 동일 DSSR 재계산 직후의 중간
> 검증 기록이다. 아래의 원고 `HOLD`는 당시 결론·문장 수정 전 상태를 뜻한다. 공개 저장소는
> signed-six 확장, 방향 감사, 변수·기준패널 민감도 분석과 결론 축소를 반영한 뒤 정리되었으며,
> 현재 결론과 재현 상태는 루트 `README.md`와 `results/generated/key_results.json`을 기준으로 한다.

- 실행일: 2026-08-04
- 검증 범위: 문제 4 — 모든 분석대상의 원자 좌표에서 동일한 프로그램으로 염기쌍 내부 변수 재계산
- 기술 상태: **PASS — 입력·실행·추출·재실행 검증 통과**
- 결론 상태: **REDUCED — ‘전반적 추가 이탈’이 아니라 ‘stretch 중심 추가 차이’만 허용**
- 원고 제출 상태: **SUPERSEDED — 2026-08-04 당시 HOLD 기록**

## 1. 이번 검증이 답한 질문

기존 원고는 여러 구조의 염기쌍 내부 변수를 좌표에서 같은 프로그램으로 다시 산출하지 않고, 보관된 값 또는 구조파일의 파생 범주를 함께 사용했다. 이번 검증에서는 공식 RCSB mmCIF 좌표를 새로 내려받아 DSSR v2.9.1 하나로 정상 기준구조와 8-oxoG 표적구조를 모두 계산했다.

이 검증은 **REFERENCE_AWARE_REPRODUCTION**이다. 기존 질문·대상 구조·예상 결과를 알고 수행했으므로 블라인드 검증이나 사전등록 연구라고 표현하지 않는다.

## 2. 입력과 도구

### 입력

- 정상 B-DNA 기준: 18개 PDB 구조
- 표적: 111D, 178D, 183D
- 공식 좌표: ASU 21개와 183D biological assembly 1개, 총 22개 파일
- 모든 입력 URL, 수집시각, ETag, 크기, SHA-256은 `01_raw_mmcif/full_v1_2026-08-04/input_manifest_full_v1.csv`에 기록했다.
- 입력 QC 결과: `PASS_INPUT_QC`

183D의 8-oxoG:C 염기쌍은 ASU 파일 안에서 완결되지 않고 결정학적 대칭 짝을 필요로 한다. 따라서 분석용 염기쌍은 RCSB assembly1에서 추출했고, 대칭으로 중복 표현되는 두 행은 모두 감사표에 남기되 한 행만 대표값으로 사용했다.

### 도구와 실행조건

- DSSR: v2.9.1-2026jul09
- 실행 파일 SHA-256: `2fbf5dd32df8a66753486b00fdb2388e4188a81f9ed7251378b2fa889b82a7dc`
- 옵션: `--more --json --nt-mapping=8OG:g`
- 분석계획 고정본 SHA-256: `e857e2cd7c4e45e0406b2f315cccd8362d5b664cdf49f89d2f8e30d915e91ece`
- 계획 수정본 SHA-256: `6e02206b9d85b862942b302af2fae2773fadabbd83349bd0af1ca012eacedd6a`
- 메타데이터 정정본 SHA-256: `f8d70ce1fbdf29fa3a86e7a2969dda2f2bb111907ab2493db3e75bd3a052e25c`

DSSR JSON의 원래 `pairs[].bp_params`를 사용했다. 여섯 값의 순서는 shear, stretch, stagger, buckle, propeller, opening이다. 염기쌍 방향을 A:T, G:C, G:A, 8OG:A, 8OG:C의 지정 방향으로 통일한 뒤 분석했다. `bp_simpleParams`는 사용하지 않았다.

## 3. 실행·재현성 검증

| 검사 | 결과 |
|---|---:|
| DSSR 입력 파일 | 22 |
| production 종료코드 0 | 22/22 |
| 동일 명령 재실행 종료코드 0 | 22/22 |
| 경로·시각 메타데이터를 제외한 수치 payload 일치 | 22/22 |
| DSSR 전체 염기쌍 행 | 264 |
| 정상 기준구조 | 18 |
| 정상 canonical Watson–Crick 염기쌍 | 230 (A:T 105, G:C 125) |
| 표적 주 분석값 | 5 |
| 표적 원시 감사값 | 6 |

첫 재현성 대조에서는 production과 rerun의 출력 경로·시작시각·종료시각까지 해시에 포함해 거짓 불일치가 발생했다. 이 결과를 삭제하지 않고 `INITIAL_FALSE_METADATA_MISMATCH` 파일로 보존한 뒤, 비수치 실행 메타데이터만 제외하여 22개 모두의 분석 수치가 동일함을 확인했다.

## 4. 좌표 재계산 결과

주 지표는 기존 질문과 맞춘 `|stretch|`, `|opening|`, `|propeller|` 3변수의 정상 기준공간 대각선 표준화 거리 `D_diagonal`이다.

| 표적 | 새 D_diagonal | 이전 D_diagonal | 새 경험 백분위 |
|---|---:|---:|---:|
| 111D G:A site 4 | 42.7414 | 44.2732 | 100.0 |
| 111D G:A site 9 | 39.8859 | 41.3253 | 100.0 |
| 178D 8-oxoG:A site 4 | 45.6972 | 47.5101 | 100.0 |
| 178D 8-oxoG:A site 9 | 41.9963 | 43.5948 | 100.0 |
| 183D 8-oxoG:C | 0.9177 | 0.9029 | 13.6 |

표적 염기쌍의 DSSR 원시 변수는 보관값과 표시 정밀도 0.001 이내에서 재현됐다. 거리의 변화는 정상 기준집합에서 발생했다. 현재 DSSR은 3IXN에서 기존 파생 범주보다 중앙 A:T 염기쌍 하나를 더 인식했고, 나머지 9개 3IXN 염기쌍에서도 과거 보관값과 차이가 있었다. 고정된 포함 규칙에 따라 현재 좌표에서 산출한 값을 주 분석에 유지했고, 3IXN 전체를 제외한 민감도 분석에서도 178D−111D 방향은 유지됐다.

## 5. 111D와 178D의 대응 자리 비교

| 자리 | 178D−111D D_diagonal | 111D 대비 반지름 거리 증가 | 두 표적 사이 직접 표준화 거리 |
|---|---:|---:|---:|
| site 4 | +2.9558 | +6.9155% | 7.7034 |
| site 9 | +2.1104 | +5.2911% | 3.9996 |

기존의 7.311%, 5.492%는 각각 6.915%, 5.291%로 갱신해야 한다. 이 비율은 **선택한 정상 기준점으로부터의 반지름 거리 증가율**이다. 두 구조의 전체 차이가 5–7%라는 뜻도 아니고 생물학적 효과크기도 아니다. 두 표적 벡터 사이의 직접 표준화 거리는 각각 7.703, 4.000이므로, 반지름 비율을 구조 전체의 유사도처럼 해석하면 안 된다.

각 자리의 원시 3변수 차이 `178D−111D`는 다음과 같다.

- site 4: `|stretch| +1.017 Å`, `|opening| −7.641°`, `|propeller| +3.836°`
- site 9: `|stretch| +0.559 Å`, `|opening| −1.943°`, `|propeller| −2.087°`

## 6. 민감도 분석과 결론 게이트

### 유지된 결과

- signed 6변수 분석: 두 자리 모두 178D−111D 거리 차이가 양수였다.
- absolute 6변수 분석: 두 자리 모두 양수였다.
- 구조 동일가중 분석: site 4 `+2.8131`, site 9 `+1.9578`이었다.
- 정상 구조 하나씩 제외: 18개 제외 경우 모두 두 자리의 방향이 양수였다.
- 3IXN 제외: site 4 `+2.9707`, site 9 `+2.1811`이었다.
- terminal pair 제외와 고해상도 구조만 사용한 분석도 방향이 양수였다.

### 결론을 제한한 결과

변수 하나씩 제외한 결과는 다음과 같다.

| 제외 변수 | site 4의 178D−111D | site 9의 178D−111D |
|---|---:|---:|
| `|stretch|` | **−2.8386** | **−0.7331** |
| `|opening|` | +7.1271 | +3.9016 |
| `|propeller|` | +2.9487 | +2.1171 |

`|stretch|`를 제외하면 두 자리 모두 방향이 뒤집혔다. 따라서 다음을 구분해야 한다.

- **허용:** 178D 8-oxoG:A에는 111D G:A와 비교해 `stretch`가 중심이 된 추가 차이가 있다.
- **불허:** 178D는 111D보다 변수 전반에서 일관되게 더 크게 이탈한다.
- **불허:** 5.3–6.9%가 두 구조 전체의 차이 또는 산화의 생물학적 효과크기다.

부트스트랩에서도 site 4의 robust 거리 차이는 95% 구간이 `−0.0240`에서 `5.5502`로 0을 아주 조금 포함했다. 주 대각선 거리의 방향은 안정적이지만 모든 거리 정의에서 완전히 확정적이라고 쓰면 안 된다.

## 7. 권장 원고 문장

> 선택한 정상 B-DNA 기준공간에서 111D의 G:A와 178D의 8-oxoG:A는 모두 극단적인 pair-internal 구조 이탈을 보였다. 따라서 큰 이탈 자체를 산화 특이적 현상으로 해석하기는 어렵다. 178D의 기준점으로부터의 반지름 거리는 대응하는 111D보다 site 4에서 6.9%, site 9에서 5.3% 컸으나, 이 방향은 stretch를 제외하면 역전되었다. 그러므로 관찰된 추가 차이는 변수 전반의 일반적 증가라기보다 stretch 중심의 차이로 해석해야 한다.

183D에 대해서는 다음 수준이 안전하다.

> 183D의 8-oxoG:C는 선택한 정상 G:C 기준공간에서 `D_diagonal=0.918`, 경험 백분위 13.6으로 기준분포 가까이에 위치했다.

## 8. 해결된 한계와 남은 한계

### 이번 실행으로 해결

- 동일 프로그램·동일 옵션으로 모든 구조의 여섯 염기쌍 내부 변수 재계산
- 원본 좌표·도구·명령·중간 JSON·파싱결과·해시 보존
- 183D 결정학적 대칭 짝의 명시적 처리
- 실행 반복에 따른 수치 재현성 확인
- 보관값과 현재 좌표 재계산값의 차이 추적

### 여전히 남음

- 111D와 178D는 각각 단일 결정구조이며, 두 자리는 독립 반복이 아니다.
- 정상 기준공간의 구조 선택 독립성과 모집단 대표성은 제한적이다.
- 구조적 차이를 생물학적 기능·돌연변이율·복구효소 인식으로 직접 연결하지 못한다.
- 관찰은 기존 공개 구조의 재해석이며 새로운 구조 측정이 아니다.
- 이 실행은 과거 v5의 208개 base-pair-step `D_step` 분석을 재현하지 않았다.

## 9. 검증 중 발견·처리한 오류

| 항목 | 심각도 | 처리 |
|---|---|---|
| 111D site 4의 최초 계획 파일 사슬 표기 오류 | 높음 | production 전 실제 DSSR endpoint로 수정하고 수정본 해시 고정 |
| 183D ASU에서 표적 쌍이 완결되지 않음 | 높음 | RCSB assembly1 사용, 대칭 중복값 모두 감사 후 대표값 1개 사용 |
| DSSR JSON의 CP1252 문자 때문에 최초 UTF-8 판독 실패 | 중간 | 실패 산출물 보존, UTF-8→CP1252 fallback을 명시적으로 적용 |
| 경로·시간 메타데이터를 포함한 최초 재실행 해시 불일치 | 중간 | 거짓 불일치 보존, 수치 payload만 비교해 22/22 일치 확인 |
| 3IXN의 과거 파생값과 현재 DSSR 값 차이 | 중간 | 현재 좌표 유래 값을 고정 규칙대로 유지하고 3IXN 제외 민감도 분석 제시 |
| 계획 파일 시각을 미래로 반올림해 기록 | 낮음 | 원본 수정 없이 별도 메타데이터 정정본 생성 |

## 10. 최종 판단

문제 4의 기술적 요구는 충족했다. 즉, **원자 좌표에서 동일 DSSR로 다시 계산하지 않았다는 한계는 이번 파이프라인 범위에서 해결됐다.** 그러나 계산이 통과했다는 사실과 원고 결론이 그대로 유지된다는 것은 다르다. 이번 결과가 허용하는 핵심 주장은 다음 두 문장이다.

1. G:A와 8-oxoG:A의 큰 pair-internal 이탈은 공통적으로 관찰되어 큰 이탈 자체는 산화 특이적이지 않다.
2. 8-oxoG:A에는 G:A보다 stretch가 중심이 된 추가 차이가 관찰되지만, 이를 변수 전반의 일반적 추가 이탈로 확대하면 안 된다.

외부 제출 전에는 학생이 같은 실행을 직접 재현하고, 표적 매핑·염기쌍 방향 변환·3IXN 차이·183D 대칭 처리·변수 제외 결과를 자신의 말로 설명할 수 있어야 한다.

## 11. 근거 파일

- 분석계획: `00_protocol/analysis_plan_lock_v1.yaml`
- 계획 수정: `00_protocol/analysis_plan_amendment_v1_1.yaml`
- 입력 manifest: `01_raw_mmcif/full_v1_2026-08-04/input_manifest_full_v1.csv`
- DSSR 실행 manifest: `02_tool_manifest/dssr_run_manifest_full_v1.csv`
- 전체 파싱값: `04_parsed_pairs/all_pairs_oriented_full_v1.csv`
- 표적값: `04_parsed_pairs/target_pairs_full_v1.csv`
- 재실행 요약: `08_logs/dssr_full_v1/dssr_run_summary_full_v1.json`
- 표적 거리: `06_analysis/target_metrics_dssr_full_v1.csv`
- 대응 자리 비교: `06_analysis/matched_comparison_dssr_full_v1.csv`
- 구조 하나씩 제외: `06_analysis/leave_one_structure_out_dssr_full_v1.csv`
- 변수 하나씩 제외: `06_analysis/leave_one_variable_out_dssr_full_v1.csv`
- 최종 결론 게이트: `06_analysis/analysis_summary_dssr_full_v1.json`
- 실행 완료 노트북: `06_analysis/coordinate_recalculation_v1.ipynb`

## 12. 방법 정의의 공식 참고문헌

- X3DNA/DSSR base-pair orientation specification: https://x3dna.org/articles/specification-of-base-pairs-in-3dna
- X3DNA simple/original base-pair parameter distinction: https://x3dna.org/highlights/details-on-the-simple-base-pair-parameters
- RCSB PDB file download service: https://www.rcsb.org/docs/programmatic-access/file-download-services
