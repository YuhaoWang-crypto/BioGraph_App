import streamlit as st
import pandas as pd
import requests
import os
import gdown

# === 页面配置 ===
st.set_page_config(
    page_title="BioGraph 全景分析台",
    page_icon="🧬",
    layout="wide"
)

# === 配置路径与ID ===
# 1. Google Drive 上的 master_graph_data.csv
DRIVE_FILE_ID = '1tM-n8EOVCPNLg9j9JfoIgg-SFYSML5XM' # 您提供的 File ID
LOCAL_CORE_PATH = "data/master_graph_data.csv"

# 2. GitHub 本地仓库里的 Swmed 文件
LOCAL_SWMED_PATH = "swmed_contacts_final_clean.csv"

# === 数据加载函数 (带缓存) ===
@st.cache_data(show_spinner=False)
def load_all_data():
    """
    加载核心数据：
    1. Core Data: 从 Google Drive 下载
    2. Swmed Data: 从本地读取
    """
    
    # --- A. 下载 Core Data ---
    os.makedirs(os.path.dirname(LOCAL_CORE_PATH), exist_ok=True)
    
    if not os.path.exists(LOCAL_CORE_PATH) or os.path.getsize(LOCAL_CORE_PATH) < 1000:
        url = f'https://drive.google.com/uc?id={DRIVE_FILE_ID}'
        # quiet=False 显示下载进度到后台日志
        gdown.download(url, LOCAL_CORE_PATH, quiet=False, fuzzy=True)
    
    try:
        df_core = pd.read_csv(LOCAL_CORE_PATH)
        # 统一列名: Source/Target -> Gene_A/Gene_B
        rename_map = {'Source': 'Gene_A', 'Target': 'Gene_B'}
        df_core.rename(columns=rename_map, inplace=True)
    except Exception as e:
        st.error(f"无法读取 master_graph_data.csv: {e}")
        return None, None, []

    # --- B. 读取 Swmed Data ---
    if os.path.exists(LOCAL_SWMED_PATH):
        try:
            df_swmed = pd.read_csv(LOCAL_SWMED_PATH)
            # 确保列名统一 (Swmed文件通常已经是 Gene_A/Gene_B 了，但也防万一)
            if 'UniProt_A' in df_swmed.columns and 'Gene_A' not in df_swmed.columns:
                 st.warning("Swmed 文件似乎缺少基因名转换，请检查文件版本。")
        except Exception as e:
            st.error(f"读取 swmed 文件出错: {e}")
            df_swmed = pd.DataFrame()
    else:
        # 如果还没上传，给个空表防止报错
        # st.warning("未找到 swmed_contacts_final_clean.csv，将只显示 BioGRID/STRING 数据。")
        df_swmed = pd.DataFrame()

    # --- C. 构建全量基因列表 ---
    genes = set()
    if not df_core.empty:
        genes.update(df_core['Gene_A'].dropna().unique())
        genes.update(df_core['Gene_B'].dropna().unique())
    
    if not df_swmed.empty:
        genes.update(df_swmed['Gene_A'].dropna().unique())
        genes.update(df_swmed['Gene_B'].dropna().unique())
        
    sorted_genes = sorted([str(g) for g in list(genes)])
    
    return df_core, df_swmed, sorted_genes

# === 辅助：UniProt 实时爬虫 ===
def fetch_uniprot_live(gene):
    try:
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            'query': f"(gene_exact:{gene}) AND (model_organism:9606) AND (reviewed:true)", 
            'fields': 'ft_region,cc_interaction', 
            'format': 'tsv', 
            'size': 1
        }
        # 伪装头
        headers = {'User-Agent': 'Mozilla/5.0'}
        r = requests.get(url, params=params, headers=headers, timeout=5)
        
        if r.status_code == 200:
            lines = r.text.strip().split('\n')
            if len(lines) > 1:
                parts = lines[1].split('\t')
                raw_region = parts[0] if len(parts) > 0 else ""
                
                res = []
                for part in raw_region.split('REGION'):
                    if 'note="' in part:
                        try:
                            loc = part.split(';')[0].strip().replace('..', '-')
                            desc = part.split('note="')[1].split('"')[0]
                            # 过滤关键词
                            if any(x in desc.lower() for x in ['interact', 'bind', 'region']):
                                res.append(f"<li><small><b>[{loc}]</b> {desc}</small></li>")
                        except: pass
                return f"<ul style='margin-bottom:0;'>{''.join(res)}</ul>" if res else "无明确位点记录"
            else:
                return "UniProt 未收录"
    except: pass
    return "网络连接失败"

# ==========================================
# 主界面逻辑
# ==========================================

st.title("🧬 BioGraph 蛋白互作全景分析台")

with st.spinner("正在初始化数据引擎..."):
    df_core, df_swmed, all_genes = load_all_data()

if df_core is None:
    st.stop() # 数据加载失败则停止

# --- 侧边栏搜索 ---
with st.sidebar:
    st.header("🔍 侦探控制台")
    default_idx = all_genes.index('TP53') if 'TP53' in all_genes else 0
    target_gene = st.selectbox(
        "选择目标蛋白", 
        options=all_genes, 
        index=default_idx,
        help="支持输入搜索"
    )
    st.markdown("---")
    st.caption(f"📚 数据库覆盖: {len(all_genes)} 个蛋白")
    st.caption("数据源: BioGRID, STRING, Swmed(AI), UniProt")

# --- 主内容区 ---
if target_gene:
    st.subheader(f"分析对象: :blue[{target_gene}]")
    
    col1, col2 = st.columns([1, 2])
    
    # 左侧：UniProt 官方信息
    with col1:
        with st.container(border=True):
            st.markdown("#### 📘 UniProt 官方位点")
            uni_html = fetch_uniprot_live(target_gene)
            st.markdown(uni_html, unsafe_allow_html=True)

    # 右侧：互作表格计算
    with col2:
        # --- 核心计算逻辑 ---
        results = {}
        
        # 1. 扫描 Core 数据 (BioGRID/STRING/AI分数)
        mask_core = (df_core['Gene_A'] == target_gene) | (df_core['Gene_B'] == target_gene)
        subset_core = df_core[mask_core]
        
        for _, row in subset_core.iterrows():
            p = row['Gene_B'] if row['Gene_A'] == target_gene else row['Gene_A']
            results[p] = {
                'Partner': p,
                'BioGRID': row.get('in_BioGRID', False) or row.get('Is_BioGRID', False), # 兼容不同列名
                'String': row.get('score_STRING', 0),
                'Swmed_Score': row.get('score_Swmed', 0),
                'My_Site': '-', 
                'Partner_Site': '-'
            }
            
        # 2. 扫描 Swmed 数据 (3D 位点)
        if not df_swmed.empty:
            mask_swm = (df_swmed['Gene_A'] == target_gene) | (df_swmed['Gene_B'] == target_gene)
            subset_swm = df_swmed[mask_swm]
            
            for _, row in subset_swm.iterrows():
                # 智能翻转逻辑
                if row['Gene_A'] == target_gene:
                    p = row['Gene_B']
                    m_site = row['Site_A_Indices']
                    p_site = row['Site_B_Indices']
                else:
                    p = row['Gene_A']
                    m_site = row['Site_B_Indices']
                    p_site = row['Site_A_Indices']
                
                # 更新或新增
                if p in results:
                    results[p]['My_Site'] = m_site
                    results[p]['Partner_Site'] = p_site
                else:
                    # Swmed 独有的新发现
                    results[p] = {
                        'Partner': p,
                        'BioGRID': False, 'String': 0, 'Swmed_Score': 0.99, # 默认高分
                        'My_Site': m_site, 'Partner_Site': p_site
                    }
        
        # 3. 转换为列表并排序
        final_data = list(results.values())
        # 排序：有3D位点的排第一 > BioGRID排第二 > 分数高排第三
        final_data.sort(key=lambda x: (x['My_Site'] != '-', x['BioGRID'], x['Swmed_Score']), reverse=True)
        
        # 4. 构造展示用 DataFrame
        display_rows = []
        for item in final_data:
            # 拼接证据标签
            tags = []
            if item['BioGRID']: tags.append("BioGRID")
            if item['String'] > 0: tags.append(f"STRING({item['String']:.2f})")
            if item['My_Site'] != '-': tags.append("Swmed(3D)")
            elif item['Swmed_Score'] > 0: tags.append(f"Swmed(Pred:{item['Swmed_Score']:.2f})")
            
            display_rows.append({
                "Partner": item['Partner'],
                "Evidence": ", ".join(tags),
                "Target_Sites": item['My_Site'],
                "Partner_Sites": item['Partner_Site']
            })
            
        df_display = pd.DataFrame(display_rows)
        
        # 渲染表格
        st.markdown(f"#### 🧬 结合伴侣与位点详情 ({len(df_display)})")
        
        st.dataframe(
            df_display,
            column_config={
                "Partner": st.column_config.TextColumn("结合伴侣", width="small"),
                "Evidence": st.column_config.TextColumn("证据来源", width="medium"),
                "Target_Sites": st.column_config.TextColumn(
                    f"{target_gene} 上的位点", 
                    help="AI 计算出的 3D 接触残基编号",
                    width="large"
                ),
                "Partner_Sites": st.column_config.TextColumn(
                    "伴侣上的位点",
                    width="large"
                )
            },
            hide_index=True,
            use_container_width=True,
            height=600
        )
