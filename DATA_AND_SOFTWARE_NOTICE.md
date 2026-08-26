# Data and software notice

## RCSB PDB 좌표

원자 좌표와 구조 메타데이터는 RCSB Protein Data Bank에서 수집했다. 저장소는 원 mmCIF 파일을
재배포하지 않고, 각 파일의 PDB ID, 공식 URL, 수집시각, 크기, SHA-256을
`config/input_manifest_full_v1.csv`에 제공한다.

결정품질 감사에 사용한 wwPDB validation XML도 재배포하지 않는다. 대신 2026-08-04에
수집한 21개 파일의 공식 URL·크기·SHA-256을
`config/validation_manifest_2026-08-04.csv`에 잠갔다. 전체 파이프라인은 새로 받은
파일을 이 manifest와 강제 대조하며, 내용이 달라지면 조용히 결과를 갱신하지 않고 중단한다.

RCSB PDB archive의 데이터는 RCSB 사용 정책에 따라 CC0로 제공된다. 다만 과학적 사용에서는
PDB와 각 원 구조 논문을 인용하는 것이 권장된다.

- RCSB usage policy: https://www.rcsb.org/pages/usage-policy
- RCSB file download services: https://www.rcsb.org/docs/programmatic-access/file-download-services

## DSSR

염기쌍 인식과 pair-internal parameter 계산에는 X3DNA-DSSR을 사용했다. DSSR은 별도 조건으로
제공되는 외부 프로그램이며 실행파일, 매뉴얼, 라이선스 파일을 이 저장소에 포함하지 않는다.
사용자는 정식 경로로 프로그램을 받아 자신의 이용조건을 확인해야 한다.

이 저장소에는 직접 작성한 호출·파싱·감사 코드와 DSSR에서 파생한 pair-level 표만 포함한다.
원 DSSR JSON은 전체 재현 실행 중 로컬에 생성되지만 기본적으로 커밋되지 않는다.

## 저장소 코드의 라이선스 상태

현재 업로드 준비본에는 코드 라이선스 파일을 의도적으로 넣지 않았다. 저장소 소유자가 공개 전
재사용 범위를 결정한 뒤 라이선스를 선택해야 한다. 라이선스를 선택하지 않으면 일반적으로 저작권
보유자가 명시적 재사용 권리를 부여하지 않은 상태가 된다.

- GitHub licensing guidance:
  https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository

## 논문과 제3자 자료

논문 PDF, 제3자 논문 원문, 학교·대회 문서, DSSR 배포물은 코드 저장소의 범위에서 제외했다.
그림은 이 저장소의 파생 데이터와 `scripts/build_dna_figures_v3.py`로 생성한 연구 결과물이다.
