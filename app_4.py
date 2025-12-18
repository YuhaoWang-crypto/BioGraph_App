import streamlit as st
import pandas as pd
import gseapy as gp
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

# --- 页面配置 ---
st.set_page_config(page_title="PPI & Drug Discovery Portal", layout="wide")

st.title("🧬 PPI & 小分子抑制剂实验设计平台")
st.markdown("支持 **通路驱动分析** 与 **自定义靶点搜索** 双模式。")

# --- 核心功能函数 ---

@st.cache_data
def get_libraries():
    """获取全量通路库名称"""
    try:
        all_libs = gp.get_library_name()
        return [l for l in all_libs if any(x in l for x in ['KEGG', 'Reactome', 'GO_Biological_Process'])]
    except:
        return ["KEGG_2021_Human", "Reactome_2022_Human"]

def get_ppi_partners(gene_symbol, limit=5):
    """从 STRING API 获取 PPI 伙伴"""
    url = "https://string-db.org/api/json/network"
    params = {"identifiers": gene_symbol, "species": 9606, "limit": limit}
    try:
        res = requests.get(url, params=params)
        data = res.json()
        partners = []
        for item in data:
            # 过滤掉查询蛋白本身
            p_name = item['preferredName_B'] if item['preferredName_A'].upper() == gene_symbol.upper() else item['preferredName_A']
            partners.append(p_name)
        return list(set(partners))
    except:
        return []

def get_inhibitors_rest(target_gene, pchembl_min=6.0):
    """使用 REST API 抓取 ChEMBL 数据"""
    search_url = f"https://www.ebi.ac.uk/chembl/api/data/target/search.json"
    search_params = {"q": target_gene}
    try:
        res = requests.get(search_url, params=search_params, timeout=10)
        targets = res.json().get('targets', [])
        
        target_id = None
        for t in targets:
            if t.get('organism') == "Homo sapiens" and t.get('target_type') == "SINGLE PROTEIN":
                target_id = t.get('target_chembl_id')
                break
        
        if not target_id: return pd.DataFrame()
        
        activity_url = f"https://www.ebi.ac.uk/chembl/api/data/activity.json"
        activity_params = {
            "target_chembl_id": target_id,
            "pchembl_value__gte": pchembl_min,
            "standard_type": "IC50",
            "order_by": "-pchembl_value",
            "limit": 10
        }
        res_act = requests.get(activity_url, params=activity_params, timeout=10)
        activities = res_act.json().get('activities', [])
        
        if not activities: return pd.DataFrame()
        
        df = pd.DataFrame(activities)
        df['target_gene'] = target_gene
        # 统一列名展示
        return df[['target_gene', 'molecule_chembl_id', 'pchembl_value', 'canonical_smiles']]
    except:
        return pd.DataFrame()

# --- 侧边栏：第一步 通路探索 ---
st.sidebar.header("第一步：通路探索")
libs = get_libraries()
selected_lib = st.sidebar.selectbox("选择通路库", libs, index=0)

@st.cache_data
def load_pathway_data(lib_name):
    return gp.get_library(lib_name)

pathway_dict = load_pathway_data(selected_lib)
pathway_list = sorted(list(pathway_dict.keys()))
selected_pathway = st.sidebar.selectbox("选择具体通路", pathway_list)

pathway_genes = []
if selected_pathway:
    pathway_genes = sorted(pathway_dict[selected_pathway])
    st.sidebar.success(f"已发现 {len(pathway_genes)} 个成员基因")
    
    with st.expander("查看当前通路基因清单"):
        st.write(", ".join(pathway_genes))

# --- 主界面：第二步 深度分析 ---
st.divider()
st.subheader("第二步：靶点/PPI 深度实验分析")

# 创建选择模式的容器
analysis_mode = st.radio(
    "选择分析模式:",
    ["从当前通路中挑选靶点", "自定义输入靶点 (直接搜索)"],
    horizontal=True
)

col1, col2 = st.columns([1, 1])

with col1:
    if analysis_mode == "从当前通路中挑选靶点":
        if not pathway_genes:
            st.warning("请先在左侧选择一个通路。")
            target_to_analyze = None
        else:
            target_to_analyze = st.selectbox("选择目标蛋白:", pathway_genes)
    else:
        target_to_analyze = st.text_input("手动输入蛋白 Symbol (如: EGFR, TP53, BRD4):", value="").upper().strip()

with col2:
    ppi_limit = st.slider("互作蛋白抓取数量", 3, 15, 5)

run_btn = st.button("🚀 开始深度分析", type="primary", use_container_width=True)

if run_btn:
    if not target_to_analyze:
        st.error("请输入或选择一个蛋白名。")
    else:
        with st.spinner(f"正在分析 {target_to_analyze} 及其 PPI 网络..."):
            # 1. 获取 PPI 伙伴
            partners = get_ppi_partners(target_to_analyze, limit=ppi_limit)
            all_targets = [target_to_analyze] + partners
            
            st.info(f"🔗 **PPI 网络节点:** {' → '.join(all_targets)}")
            
            # 2. 批量获取抑制剂
            all_inhibitor_data = []
            progress_bar = st.progress(0)
            
            for i, t in enumerate(all_targets):
                df_inh = get_inhibitors_rest(t)
                if not df_inh.empty:
                    all_inhibitor_data.append(df_inh)
                progress_bar.progress((i + 1) / len(all_targets))
            
            if all_inhibitor_data:
                final_df = pd.concat(all_inhibitor_data, ignore_index=True)
                final_df['pchembl_value'] = pd.to_numeric(final_df['pchembl_value'], errors='coerce')
                final_df = final_df.sort_values('pchembl_value', ascending=False).dropna(subset=['pchembl_value'])
                
                # --- 展示表格 ---
                st.write("### 💊 综合抑制剂清单 (目标及其互作蛋白)")
                st.dataframe(final_df, use_container_width=True)
                
                # --- 展示结构图 (Top 3 活性分子) ---
                st.write("### 🎨 高活性分子结构参考 (Top 3)")
                top_3 = final_df.head(3)
                img_cols = st.columns(3)
                
                for idx, row in enumerate(top_3.itertuples()):
                    with img_cols[idx]:
                        mol = Chem.MolFromSmiles(row.canonical_smiles)
                        if mol:
                            # 增加绘制细节
                            img = Draw.MolToImage(mol, size=(300, 300))
                            st.image(img, caption=f"{row.target_gene}: {row.molecule_chembl_id}\n(pChEMBL: {row.pchembl_value})")
                
                # 下载
                csv = final_df.to_csv(index=False).encode('utf-8')
                st.download_button(
                    label="📥 下载分析结果 (CSV)",
                    data=csv,
                    file_name=f"{target_to_analyze}_PPI_Drug_Analysis.csv",
                    mime="text/csv",
                )
            else:
                st.warning(f"在 ChEMBL 中未找到 {target_to_analyze} 及其伙伴的高活性抑制剂数据。")

# 底部页脚
st.markdown("---")
st.caption("数据来源: STRING (PPI) | ChEMBL (小分子) | Enrichr (通路库) | RDKit (结构绘制)")
