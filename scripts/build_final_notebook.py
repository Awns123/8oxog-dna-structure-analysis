"""Build the reader-facing final signed-six reproduction notebook."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "notebooks" / "final_signed_six_reproduction.ipynb"


def markdown(source: str) -> dict[str, object]:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(keepends=True)}


def code(source: str) -> dict[str, object]:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source.splitlines(keepends=True),
    }


def main() -> int:
    cells = [
        markdown(
            """# 8-oxoG:A matched-structure analysis

## tl;dr

111D의 G:A와 178D의 8-oxoG:A를 동일한 signed-six 기준공간에서 비교한다. 두 구조 모두 정상 A:T
기준에서 크게 이탈하지만, 178D의 추가 radial difference는 stretch에 강하게 의존하며 stretch
제외 시 두 위치 모두 방향이 역전된다.
"""
        ),
        markdown(
            """## Context & Methods

여섯 pair-internal rigid-body parameter(shear, stretch, stagger, buckle, propeller, opening)를 염기 순서에
맞게 방향 정규화한 뒤, 정상 B-DNA 기준패널의 평균과 표본 표준편차로 표준화한다. 각 표적의 radial
distance와 111D–178D 사이의 direct standardized distance를 구분한다.

### Key Assumptions

- 111D와 178D는 각각 하나의 결정구조이며 site 4·9는 독립 반복이 아니다.
- 거리값은 선택된 기준패널과 변수 정의에 조건부이다.
- 분석은 구조 사례의 재해석이며 산화의 모집단 인과효과를 추정하지 않는다.
"""
        ),
        code(
            """from pathlib import Path
import sys

repo = Path.cwd()
if not (repo / 'run_analysis.py').is_file():
    repo = repo.parent
sys.path.insert(0, str(repo))

from run_analysis import analyze

result = analyze(output_dir=None)
summary = result['summary']
values = summary['values']
print(summary['status'])
"""
        ),
        markdown(
            """## Data

입력은 DSSR 결과에서 추출하고 방향 정규화한 pair-level CSV이다. 정상 기준패널은 A:T 105쌍과
G:C 125쌍이며, 표적은 111D·178D의 두 대응 위치와 별도 183D 8-oxoG:C 사례이다.
"""
        ),
        code(
            """summary['reference_pairs'], summary['independent_target_structures']
"""
        ),
        markdown("## Results\n\n### Radial and direct comparisons\n"),
        code("result['matched_comparison'].round(4)\n"),
        markdown("### Stretch dependence\n"),
        code(
            """print(f"site 4: stretch direct D^2 share = {values['site4_stretch_direct_share_pct']:.2f}%")
print(f"site 9: stretch direct D^2 share = {values['site9_stretch_direct_share_pct']:.2f}%")
print(f"omit stretch delta D: site 4 = {values['site4_omit_stretch_delta_D']:+.3f}")
print(f"omit stretch delta D: site 9 = {values['site9_omit_stretch_delta_D']:+.3f}")
"""
        ),
        markdown(
            """## Takeaways

1. 큰 pair-internal 이탈의 존재 자체는 8-oxoG:A에 고유하지 않다.
2. 이 matched case에서 178D는 111D보다 radial distance가 더 크지만, 그 방향은 stretch 포함 여부에 민감하다.
3. 따라서 허용되는 결론은 전반적 산화 효과가 아니라 **stretch 중심의 조건부 추가 차이**이다.

![signed-six reference distance](../results/figures/figure1_signed6_reference_distance_v3.png)
"""
        ),
    ]
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(notebook, ensure_ascii=False, indent=1) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(f"WROTE {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
