# synplan_route_cluster_utils.py
"""
SynPlanner route CGR 생성 + 클러스터링 헬퍼 모듈

이 파일 하나만 import 하면, 아래 코드들이 그대로 동작하도록 정리되어 있음.

- compose_all_route_cgrs
- reduce_cgr_to_main_component
- compose_all_reduced_route_cgrs
- cluster_routes
- subcluster_all_clusters
- post_process_subgroup
- cgr_display
- get_route_svg_from_json
- routes_clustering_report
- routes_subclustering_report
- display, HTML, SVG (IPython.display)

Colab 예시:
    from synplan_route_cluster_utils import *

    # all_route_cgrs = compose_all_route_cgrs(routes_dict_2)
    # all_reduced_route_cgrs = compose_all_reduced_route_cgrs(all_route_cgrs)
    # clusters = cluster_routes(all_reduced_route_cgrs, use_strat=False)
    # ...
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

# ========== SynPlanner import들 ==========
from synplan.chem.reaction_routes.route_cgr import compose_all_route_cgrs
from synplan.chem.reaction_routes.clustering import (
    cluster_routes,
    subcluster_all_clusters,
    post_process_subgroup,
)
from synplan.chem.reaction_routes.visualisation import cgr_display
from synplan.utils.visualisation import (
    get_route_svg_from_json,
    routes_clustering_report,
    routes_subclustering_report,
)

# ========== Jupyter / Colab 시각화 ==========
from IPython.display import display, HTML, SVG


# ========== Reduced CGR 유틸 함수들 (노트북에서 쓰던 것 그대로) ==========

def reduce_cgr_to_main_component(route_cgr):
    """
    하나의 Route CGR에서 가장 큰 connected component만 남긴
    'reduced' CGR을 생성한다.

    Parameters
    ----------
    route_cgr : CGR 객체 (예: CGRTools CondensedGraphOfReaction)
        compose_route_cgr / compose_all_route_cgrs 의 value로 나오는 객체.

    Returns
    -------
    reduced_cgr : CGR 객체
        가장 큰 component만 남긴 새로운 CGR.
        connected_components 가 비어 있으면 원본 CGR을 그대로 반환.
    """
    components = getattr(route_cgr, "connected_components", None)
    if not components:
        # component 정보가 없거나 비어 있으면 그냥 원본 반환
        return route_cgr

    # 각 component 는 원자 인덱스 집합이라고 가정하고 크기가 가장 큰 것 선택
    main_component = max(components, key=len)

    # 해당 component만 남긴 sub-CGR 생성
    reduced_cgr = route_cgr.substructure(main_component)
    return reduced_cgr


def compose_all_reduced_route_cgrs(tree_or_cgrs, route_id: Optional[str] = None) -> Dict[str, Any]:
    """
    여러 Route CGR에 대해 'reduced CGR'를 일괄 생성한다.
    - dict 기반({route_id: CGR} 또는 Tree 기반 둘 다 지원.

    Parameters
    ----------
    tree_or_cgrs : dict 또는 synplan.mcts.tree.Tree
        1) dict 인 경우: {route_id: CGR} 형태 (compose_all_route_cgrs 결과)
        2) Tree 인 경우: winning routes 를 내부에서 CGR로 변환 후 사용.
    route_id : optional
        특정 route 하나만 대상으로 하고 싶을 때 route_id 지정.
        None 이면 모든 route 에 대해 reduced CGR 생성.

    Returns
    -------
    dict
        {route_id: reduced_CGR} 형태의 딕셔너리.
    """
    # 1) 이미 {route_id: CGR} dict 를 받은 경우
    if isinstance(tree_or_cgrs, dict):
        route_cgrs = tree_or_cgrs
    else:
        # 2) Tree 를 받은 경우: 먼저 모든 route 에 대해 CGR 생성
        route_cgrs = compose_all_route_cgrs(tree_or_cgrs, route_id=route_id)

    def _single(rid: str):
        cgr = route_cgrs.get(rid)
        if cgr is None:
            return None
        return reduce_cgr_to_main_component(cgr)

    # 특정 route_id 만 요청한 경우
    if route_id is not None:
        if route_id not in route_cgrs:
            raise KeyError(f"Route {route_id} not in provided CGR dict.")
        reduced = _single(route_id)
        return {route_id: reduced} if reduced is not None else {}

    # 전체 route 에 대해 reduced CGR 생성
    reduced_dict: Dict[str, Any] = {}
    for rid in sorted(route_cgrs):
        reduced = _single(rid)
        if reduced is not None:
            reduced_dict[rid] = reduced

    return reduced_dict

