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

st.title("🧬 PPI & 小分子抑制剂实验设计平台 (稳定版)")
st.markdown("该版本直接连接 ChEMBL REST API，避开了客户端连接错误。")

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
            p_name = item['preferredName_B'] if item['preferredName_A'].upper() == gene_symbol.upper() else item['preferredName_A']
            partners.append(p_name)
        return list(set(partners))
    except:
        return []

def get_inhibitors_rest(target_gene, pchembl_min=6.0):
    """使用 REST API 直接抓取 ChEMBL 抑制剂数据 (更稳定)"""
    # 1. 搜索靶点
    search_url = f"https://www.ebi.ac.uk/chembl/api/data/target/search.json"
    search_params = {"q": target_gene}
    try:
        res = requests.get(search_url, params=search_params, timeout=10)
        targets = res.json().get('targets', [])
        
        # 筛选人类单蛋白靶点
        target_id = None
        for t in targets:
            if t.get('organism') == "Homo sapiens" and t.get('target_type') == "SINGLE PROTEIN":
                target_id = t.get('target_chembl_id')
                break
        
        if not target_id: return pd.DataFrame()
        
        # 2. 抓取活性数据
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
        return df[['target_gene', 'molecule_chembl_id', 'pchembl_value', 'canonical_smiles']]
    except Exception as e:
        # st.error(f"API Error for {target_gene}: {e}")
        return pd.DataFrame()

# --- 侧边栏：通路探索 ---
st.sidebar.header("第一步：通路搜索")
libs = get_libraries()
selected_lib = st.sidebar.selectbox("选择通路库", libs, index=0)

@st.cache_data
def load_pathway_data(lib_name):
    return gp.get_library(lib_name)

pathway_dict = load_pathway_data(selected_lib)
pathway_list = sorted(list(pathway_dict.keys()))
selected_pathway = st.sidebar.selectbox("选择具体通路", pathway_list)

if selected_pathway:
    genes = sorted(pathway_dict[selected_pathway])
    st.sidebar.success(f"发现 {len(genes)} 个成员基因")
    
    # --- 主界面 ---
    st.subheader(f"📍 通路成员解析: {selected_pathway}")
    st.text_area("基因清单 (可复制)", ", ".join(genes), height=100)
    
    st.divider()
    st.subheader("第二步：靶点深度实验分析")
    
    col1, col2 = st.columns([1, 2])
    with col1:
        target_input = st.selectbox("从通路中选择靶点", genes)
        ppi_limit = st.slider("互作蛋白抓取数量", 3, 10, 5)
        run_btn = st.button("开始深度分析", type="primary")

    if run_btn:
        with st.spinner(f"正在分析 {target_input} 及其 PPI 网络..."):
            # 1. PPI 获取
            partners = get_ppi_partners(target_input, limit=ppi_limit)
            all_targets = [target_input] + partners
            st.info(f"🔗 **PPI 网络节点:** {', '.join(all_targets)}")
            
            # 2. 抑制剂获取
            all_res = []
            for t in all_targets:
                df_inh = get_inhibitors_rest(t)
                if not df_inh.empty:
                    all_res.append(df_inh)
            
            if all_res:
                final_df = pd.concat(all_res, ignore_index=True)
                final_df['pchembl_value'] = pd.to_numeric(final_df['pchembl_value'])
                final_df = final_df.sort_values('pchembl_value', ascending=False)
                
                st.write("### 💊 综合抑制剂清单")
                st.dataframe(final_df, use_container_width=True)
                
                st.write("### 🎨 活性最强分子结构预览")
                top_smiles = final_df.iloc[0]['canonical_smiles']
                mol = Chem.MolFromSmiles(top_smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 400))
                    st.image(img, caption=f"Top Molecule: {final_df.iloc[0]['molecule_chembl_id']}")
                
                csv = final_df.to_csv(index=False).encode('utf-8')
                st.download_button("下载 CSV 结果", csv, "ppi_inhibitors.csv", "text/csv")
            else:
                st.warning("未能在线找到高活性小分子。")
