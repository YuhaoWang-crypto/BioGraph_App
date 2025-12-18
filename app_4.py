import streamlit as st
import pandas as pd
import gseapy as gp
import requests
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

# --- 页面配置 ---
st.set_page_config(page_title="PPI & Drug Discovery Portal", layout="wide")

st.title("🧬 PPI & 小分子抑制剂实验设计平台")
st.markdown("通过通路搜索靶点，实时挖掘蛋白质相互作用网络及其潜在小分子抑制剂。")

# --- 初始化 ChEMBL 客户端 ---
target_api = new_client.target
activity_api = new_client.activity

# --- 核心功能函数 ---

@st.cache_data
def get_libraries():
    """获取全量通路库名称"""
    all_libs = gp.get_library_name()
    return [l for l in all_libs if any(x in l for x in ['KEGG', 'Reactome', 'GO_Biological_Process'])]

def get_ppi_partners(gene_symbol, limit=5):
    """从 STRING API 获取 PPI 伙伴"""
    string_api_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": gene_symbol,
        "species": 9606,  # Human
        "limit": limit
    }
    try:
        res = requests.get(string_api_url, params=params)
        data = res.json()
        partners = []
        for item in data:
            p_name = item['preferredName_B'] if item['preferredName_A'].upper() == gene_symbol.upper() else item['preferredName_A']
            partners.append(p_name)
        return list(set(partners))
    except:
        return []

def get_inhibitors(target_gene, pchembl_min=6.0):
    """从 ChEMBL 获取抑制剂数据"""
    try:
        # 搜索靶点
        targets = target_api.search(target_gene).filter(organism="Homo sapiens", target_type="SINGLE PROTEIN")
        if not targets: return pd.DataFrame()
        
        target_id = targets[0]['target_chembl_id']
        
        # 抓取活性数据 (Top 10)
        activities = activity_api.filter(
            target_chembl_id=target_id, 
            pchembl_value__gte=pchembl_min,
            standard_type="IC50"
        ).order_by('-pchembl_value')[:10]
        
        df = pd.DataFrame(list(activities))
        if not df.empty:
            df['target_gene'] = target_gene
            return df[['target_gene', 'molecule_chembl_id', 'pchembl_value', 'canonical_smiles']]
    except:
        pass
    return pd.DataFrame()

# --- 侧边栏：通路探索 ---
st.sidebar.header("第一步：通路搜索")
libs = get_libraries()
selected_lib = st.sidebar.selectbox("选择通路库", libs, index=0)

# 加载选中库的所有通路
@st.cache_data
def load_pathway_data(lib_name):
    return gp.get_library(lib_name)

pathway_dict = load_pathway_data(selected_lib)
pathway_list = sorted(list(pathway_dict.keys()))

selected_pathway = st.sidebar.selectbox("选择具体通路", pathway_list)

if selected_pathway:
    genes = sorted(pathway_dict[selected_pathway])
    st.sidebar.success(f"发现 {len(genes)} 个成员基因")
    
    # --- 主界面：解析基因 ---
    st.subheader(f"📍 通路成员解析: {selected_pathway}")
    st.text_area("基因清单 (可复制)", ", ".join(genes), height=100)
    
    # --- 第二步：靶点分析 ---
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
            st.write(f"🔗 **PPI 网络节点:** {', '.join(all_targets)}")
            
            # 2. 抑制剂获取
            all_res = []
            for t in all_targets:
                df_inh = get_inhibitors(t)
                if not df_inh.empty:
                    all_res.append(df_inh)
            
            if all_res:
                final_df = pd.concat(all_res, ignore_index=True)
                final_df = final_df.sort_values('pchembl_value', ascending=False)
                
                # --- 结果显示 ---
                st.write("### 💊 综合抑制剂清单 (实时在线数据)")
                st.dataframe(final_df, use_container_width=True)
                
                # --- 结构显示 ---
                st.write("### 🎨 活性最强分子结构预览 (Top 1)")
                top_smiles = final_df.iloc[0]['canonical_smiles']
                mol = Chem.MolFromSmiles(top_smiles)
                if mol:
                    img = Draw.MolToImage(mol, size=(400, 400))
                    st.image(img, caption=f"Top Molecule: {final_df.iloc[0]['molecule_chembl_id']}")
                
                # 下载按钮
                csv = final_df.to_csv(index=False).encode('utf-8')
                st.download_button("下载结果表格 (CSV)", csv, "ppi_inhibitors.csv", "text/csv")
            else:
                st.warning("未能找到相关小分子抑制剂。")
