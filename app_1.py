import streamlit as st
import pandas as pd
import networkx as nx
import os
from pyvis.network import Network
import streamlit.components.v1 as components
import gdown
import itertools
from collections import Counter

# ==========================================
# ⚙️ 1. 页面配置与基础设置
# ==========================================
st.set_page_config(page_title="BioGraph 全景分析平台", layout="wide", page_icon="🧬")

# 路径设置
DATA_FILE = "data/master_graph_data.csv"
CACHE_DIR = "checkpoints"
os.makedirs(CACHE_DIR, exist_ok=True)

# === 硬编码的功能模块字典 (来自你的 Colab) ===
FUNCTION_THEMES = {
    '🛡️ 免疫系统 (Immune)': ['Immunoglobulin', 'Cytokine', 'MHC', 'IPR007110', 'Interferon', 'Inflammation', 'T-cell', 'B-cell'],
    '⚡ 代谢与能量 (Metabolism)': ['Mitochondrion', 'ATPase', 'Kinase', 'Glycolysis', 'Lipid', 'IPR000719', 'IPR004033'],
    '🧬 基因表达 (Gene Reg)': ['Nucleus', 'Zinc Finger', 'Homeobox', 'transcription', 'Histone', 'Chromatin', 'RNA-binding'],
    '🏗️ 细胞骨架 (Cytoskeleton)': ['Actin', 'Tubulin', 'Kinesin', 'Microtubule', 'Myosin', 'Filament', 'Centrosome'],
    '📡 信号传导 (Signaling)': ['Receptor', 'Membrane', 'SH2', 'G-protein', 'Phosphatase', 'Growth factor', 'MAPK'],
    '💀 细胞死亡 (Cell Death)': ['Apoptosis', 'Caspase', 'Bcl-2', 'Death domain', 'Autophagy', 'Lysosome', 'Necrosis'],
    '🧠 神经系统 (Neuro)': ['Synapse', 'Axon', 'Channel', 'Neurotransmitter', 'GABA', 'Glutamate'],
    '🧊 自噬与降解 (Autophagy)': ['Autophagy', 'Lysosome', 'Ubiquitin', 'Proteasome', 'LC3', 'SQSTM1', 'ATG'],
    '🔥 缺氧与应激 (Hypoxia/Stress)': ['Hypoxia', 'HIF', 'HSPA5', 'UPR', 'ER stress', 'Heat shock'],
    '🧬 DNA损伤修复 (DNA Repair)': ['BRCA', 'RAD51', 'ATM', 'ATR', 'Damage', 'Repair', 'Checkpoint']
}

# ==========================================
# 📥 2. 数据加载核心
# ==========================================
@st.cache_resource
def load_graph_data():
    # 确保 data 文件夹存在
    os.makedirs(os.path.dirname(DATA_FILE), exist_ok=True)

    # 下载逻辑
    if not os.path.exists(DATA_FILE) or os.path.getsize(DATA_FILE) < 1000:
        file_id = '1tM-n8EOVCPNLg9j9JfoIgg-SFYSML5XM'  # 你的 File ID
        url = f'https://drive.google.com/uc?id={file_id}'
        with st.spinner('正在从云端下载数据...'):
            gdown.download(url, DATA_FILE, quiet=False, fuzzy=True)

    # 读取逻辑
    with st.spinner('正在构建内存图谱...'):
        if not os.path.exists(DATA_FILE):
            return None, None
        
        try:
            df = pd.read_csv(DATA_FILE)
            # 兼容性处理：如果 CSV 里没有 anno 列，创建一个空的
            if 'Source_Anno' not in df.columns: df['Source_Anno'] = ''
            if 'Target_Anno' not in df.columns: df['Target_Anno'] = ''
            # 兼容性处理：IDR 列
            if 'Source_IDR' not in df.columns: df['Source_IDR'] = -1.0 # 假设
            
        except Exception as e:
            st.error(f"文件读取错误: {e}")
            os.remove(DATA_FILE)
            return None, None

        G = nx.Graph()
        for _, row in df.iterrows():
            s, t = row['Source'], row['Target']
            
            # 读取边属性
            score = row.get('score_Swmed', 0.0)
            is_known = row.get('in_BioGRID', False)
            
            # 写入边
            G.add_edge(s, t, swmed=score, known=is_known, resistance=1.0/max(score, 0.001))
            
            # 写入节点属性 (如果节点已存在则更新，简单起见这里覆盖)
            # 注意：CSV 是边列表，节点属性需要去重。这里为了速度，每次都写一遍也没事
            for node, prefix in [(s, 'Source'), (t, 'Target')]:
                G.add_node(node, 
                           loc=str(row.get(f'{prefix}_Loc', 'Unknown')),
                           dom=str(row.get(f'{prefix}_Domain', '')),
                           anno=str(row.get(f'{prefix}_Anno', '')),
                           # 尝试读取 IDR (如果你的 CSV 后来补充了这列)
                           idr=float(row.get(f'{prefix}_IDR', -1.0))
                )
        return G, df

G, master_df = load_graph_data()

if G is None:
    st.error("数据加载失败，请刷新页面。")
    st.stop()

ALL_GENES = sorted(list(G.nodes()))

# ==========================================
# 🎨 3. 侧边栏与模块选择
# ==========================================
st.sidebar.title("🧬 BioGraph Pro")
module = st.sidebar.radio("选择功能模块:", [
    "🕵️‍♂️ 深度蛋白侦探 (Detective)",
    "📚 GO/功能 智能搜索 (GO Search)",
    "⚖️ 全景双轨分析 (Panorama)",
    "🔀 模块串扰挖掘 (Crosstalk)",
    "🗺️ 定向路径挖掘 (Pathfinder)",
    "♟️ 战略阻断模拟 (Blockade)",
    "🧪 组学幕后黑手 (Omics Miner)",
    "💧 遗失的魔术贴 (IDR/LLPS)"
])

# ==========================================
# 🕵️‍♂️ 模块 1: 深度蛋白侦探
# ==========================================
if module.startswith("🕵️‍♂️"):
    st.header("🕵️‍♂️ 深度蛋白侦探")
    col1, col2, col3 = st.columns([2, 1, 1])
    target = col1.selectbox("目标蛋白:", ALL_GENES, index=ALL_GENES.index('TP53') if 'TP53' in ALL_GENES else 0)
    min_score = col2.slider("AI 分数阈值:", 0.0, 1.0, 0.0)
    filter_loc = col3.text_input("筛选定位 (可选):", placeholder="如 Nucleus")

    if st.button("开始侦查"):
        props = G.nodes[target]
        st.info(f"**{target}** 信息: 📍 {props.get('loc')} | 🧩 {props.get('dom')[:50]}...")
        
        data = []
        for n in G.neighbors(target):
            edge = G[target][n]
            n_props = G.nodes[n]
            
            # 筛选
            if edge.get('swmed', 0) < min_score and not edge.get('known'): continue
            if filter_loc and filter_loc.lower() not in n_props.get('loc', '').lower(): continue
            
            data.append({
                "Partner": n,
                "Type": "✅已知" if edge.get('known') else "🤖预测",
                "Score": edge.get('swmed', 0),
                "Location": n_props.get('loc', 'Unk'),
                "Annotation": n_props.get('anno', '')[:50] + "..."
            })
        
        df_res = pd.DataFrame(data).sort_values("Score", ascending=False)
        st.write(f"找到 {len(df_res)} 个伙伴:")
        st.dataframe(df_res, use_container_width=True)

# ==========================================
# 📚 模块 2: GO/功能 智能搜索
# ==========================================
elif module.startswith("📚"):
    st.header("📚 GO/功能 智能搜索")
    st.markdown("通过基因本体论 (GO) 或关键词搜索具有特定功能的蛋白集合。")
    
    col1, col2 = st.columns([3, 1])
    keywords = col1.text_input("输入关键词 (逗号分隔):", "autophagy, hypoxia")
    match_mode = col2.radio("匹配模式:", ["任意匹配 (OR)", "全部匹配 (AND)"])
    
    if st.button("全库搜索"):
        kws = [k.strip().lower() for k in keywords.split(',') if k.strip()]
        if not kws: st.warning("请输入关键词")
        else:
            hits = []
            for n, d in G.nodes(data=True):
                anno = str(d.get('anno', '')).lower()
                # 同时也搜名字
                name_match = any(kw in n.lower() for kw in kws)
                
                if match_mode.startswith("全部"):
                    is_match = all(kw in anno for kw in kws)
                else:
                    is_match = any(kw in anno for kw in kws) or name_match
                
                if is_match:
                    hits.append({
                        "Gene": n,
                        "Loc": d.get('loc'),
                        "Snippet": anno[:100] + "..."
                    })
            
            df_hits = pd.DataFrame(hits)
            st.success(f"🔍 找到 {len(df_hits)} 个相关蛋白")
            if not df_hits.empty:
                st.dataframe(df_hits, use_container_width=True)

# ==========================================
# ⚖️ 模块 3: 全景双轨分析
# ==========================================
elif module.startswith("⚖️"):
    st.header("⚖️ 全景双轨分析 (Balanced View)")
    st.markdown("对比 **已知通路 ** 与 **AI 预测 **，发现模块间的潜在联系。")
    
    opts = list(FUNCTION_THEMES.keys())
    c1, c2 = st.columns(2)
    theme_a = c1.selectbox("模块 A:", opts, index=0)
    theme_b = c2.selectbox("模块 B:", opts, index=1)
    
    if st.button("运行全景扫描"):
        if theme_a == theme_b:
            st.warning("请选择两个不同的模块。")
        else:
            # 1. 筛选节点
            def get_nodes(theme):
                kws = FUNCTION_THEMES[theme]
                res = set()
                for n, d in G.nodes(data=True):
                    content = (str(d.get('loc')) + str(d.get('dom')) + str(d.get('anno'))).lower()
                    if any(k.lower() in content for k in kws): res.add(n)
                return res
            
            nodes_a = get_nodes(theme_a)
            nodes_b = get_nodes(theme_b)
            
            st.write(f"🔹 模块 A ({theme_a}): {len(nodes_a)} 节点")
            st.write(f"🔸 模块 B ({theme_b}): {len(nodes_b)} 节点")
            
            # 2. 找连接
            edges = []
            for u in nodes_a:
                for v in nodes_b:
                    if G.has_edge(u, v):
                        e = G[u][v]
                        edges.append({
                            "Source (A)": u, "Target (B)": v,
                            "Type": "✅已知" if e.get('known') else "🚀预测",
                            "Score": e.get('swmed', 0)
                        })
            
            df = pd.DataFrame(edges).sort_values("Score", ascending=False)
            
            tab1, tab2 = st.tabs(["✅ 已知通路 (Benchmark)", "🚀 潜在发现 (Discovery)"])
            with tab1:
                st.dataframe(df[df['Type'].str.contains("已知")], use_container_width=True)
            with tab2:
                st.dataframe(df[df['Type'].str.contains("预测")], use_container_width=True)

# ==========================================
# 🔀 模块 4: 模块串扰挖掘
# ==========================================
elif module.startswith("🔀"):
    st.header("🔀 功能模块串扰挖掘 (Crosstalk)")
    st.markdown("寻找连接两个功能模块的 **关键桥梁蛋白 (Bridge Proteins)**。")
    
    opts = list(FUNCTION_THEMES.keys())
    c1, c2 = st.columns(2)
    theme_a = c1.selectbox("模块 A:", opts, index=2)
    theme_b = c2.selectbox("模块 B:", opts, index=4)
    
    if st.button("寻找桥梁"):
        # 复用筛选逻辑
        def get_nodes(theme):
            kws = FUNCTION_THEMES[theme]
            return {n for n, d in G.nodes(data=True) 
                    if any(k.lower() in (str(d.get('loc'))+str(d.get('dom'))+str(d.get('anno'))).lower() for k in kws)}
        
        na, nb = get_nodes(theme_a), get_nodes(theme_b)
        
        # 寻找既连接A又连接B的中间人
        bridges = {}
        for n in G.nodes():
            if n in na or n in nb: continue # 排除模块内部人员（寻找外部桥梁）
            
            neighbors = set(G.neighbors(n))
            links_a = len(neighbors.intersection(na))
            links_b = len(neighbors.intersection(nb))
            
            if links_a > 0 and links_b > 0:
                bridges[n] = links_a * links_b # 简单的打分：连接数的乘积
        
        if bridges:
            sorted_bridges = sorted(bridges.items(), key=lambda x: x[1], reverse=True)[:20]
            data = [{"Bridge Protein": k, "Links to A": "Many", "Links to B": "Many", "Score": v, 
                     "Loc": G.nodes[k].get('loc')} for k, v in sorted_bridges]
            st.success(f"找到 {len(bridges)} 个潜在桥梁，展示 Top 20:")
            st.dataframe(pd.DataFrame(data))
        else:
            st.warning("未找到明显的外部桥梁。")

# ==========================================
# 🗺️ 模块 5: 定向路径挖掘
# ==========================================
elif module.startswith("🗺️"):
    st.header("🗺️ 定向路径挖掘机")
    c1, c2, c3 = st.columns(3)
    start = c1.selectbox("起点:", ALL_GENES, index=0)
    end = c2.selectbox("终点:", ALL_GENES, index=1)
    via = c3.multiselect("必经点 (可选):", ALL_GENES)
    
    if st.button("规划路线"):
        try:
            path = []
            if not via:
                path = nx.shortest_path(G, start, end, weight='resistance')
                cost = nx.shortest_path_length(G, start, end, weight='resistance')
            else:
                # 简单实现：经过第一个必经点
                v = via[0]
                p1 = nx.shortest_path(G, start, v, weight='resistance')
                p2 = nx.shortest_path(G, v, end, weight='resistance')
                path = p1 + p2[1:]
                cost = "N/A"
            
            st.code(" -> ".join(path))
            st.caption(f"路径总阻力: {cost}")
            
            # 详情
            steps = []
            for i in range(len(path)-1):
                u, v = path[i], path[i+1]
                e = G[u][v]
                steps.append({"Step": f"{u}->{v}", "Type": "Known" if e.get('known') else "AI", "Score": e.get('swmed')})
            st.dataframe(pd.DataFrame(steps))
            
        except nx.NetworkXNoPath:
            st.error("无通路")

# ==========================================
# ♟️ 模块 6: 战略阻断模拟
# ==========================================
elif module.startswith("♟️"):
    st.header("♟️ 战略阻断模拟")
    c1, c2 = st.columns(2)
    s = c1.selectbox("起点:", ALL_GENES, key='bk_s')
    e = c2.selectbox("终点:", ALL_GENES, key='bk_e')
    
    if st.button("推演代偿"):
        try:
            base_path = nx.shortest_path(G, s, e, weight='resistance')
            base_cost = nx.shortest_path_length(G, s, e, weight='resistance')
            
            res = []
            targets = base_path[1:-1]
            if not targets: st.warning("路径太短，无中间靶点。")
            
            for t in targets:
                G_tmp = G.copy()
                G_tmp.remove_node(t)
                try:
                    new_cost = nx.shortest_path_length(G_tmp, s, e, weight='resistance')
                    impact = (new_cost - base_cost)/base_cost * 100
                    res.append({"Target": t, "Impact": impact, "New Cost": new_cost})
                except:
                    res.append({"Target": t, "Impact": 9999, "New Cost": "Inf"})
            
            st.dataframe(pd.DataFrame(res).sort_values("Impact", ascending=False).style.apply(
                lambda x: ['background-color: #ffcccc' if x['Impact'] > 1000 else '' for i in x], axis=1
            ))
        except:
            st.error("计算出错或无通路")

# ==========================================
# 🧪 模块 7: 组学幕后黑手
# ==========================================
elif module.startswith("🧪"):
    st.header("🧪 组学数据挖掘机 (Omics Miner)")
    st.markdown("输入一组看似不相关的差异蛋白，寻找连接它们的**隐形枢纽 (Hidden Hubs)**。")
    
    default_list = "TP53, EGFR, BCL2, CDK1, MYC"
    user_input = st.text_area("输入基因列表 (逗号或换行分隔):", default_list)
    
    if st.button("挖掘黑手"):
        genes = [g.strip().upper() for g in user_input.replace('\n', ',').split(',') if g.strip()]
        valid = [g for g in genes if g in G]
        st.write(f"✅ 有效识别: {len(valid)} 个基因")
        
        if len(valid) < 2:
            st.error("基因太少，无法分析")
        else:
            # 寻找隐形枢纽
            hubs = {}
            for n in G.nodes():
                if n in valid: continue
                neighbors = set(G.neighbors(n))
                overlap = len(neighbors.intersection(set(valid)))
                if overlap >= 2:
                    hubs[n] = overlap
            
            if hubs:
                top_hubs = sorted(hubs.items(), key=lambda x: x[1], reverse=True)[:10]
                df_hubs = pd.DataFrame(top_hubs, columns=["Hidden Hub", "Connected Targets"])
                st.success("🕵️‍♂️ 发现潜在的幕后调控者:")
                st.dataframe(df_hubs)
                
                # 可视化 Top 1
                top_hub = top_hubs[0][0]
                sub_nodes = valid + [top_hub]
                subG = G.subgraph(sub_nodes)
                
                net = Network(height="500px", width="100%", bgcolor="#222")
                for n in subG.nodes():
                    color = "#ff0000" if n == top_hub else "#00ff00"
                    net.add_node(n, label=n, color=color)
                for u, v in subG.edges():
                    net.add_edge(u, v, color="gray")
                
                path = "tmp_omics.html"
                net.save_graph(path)
                with open(path, 'r', encoding='utf-8') as f:
                    components.html(f.read(), height=520)

# ==========================================
# 💧 模块 8: 遗失的魔术贴 (IDR/LLPS)
# ==========================================
elif module.startswith("💧"):
    st.header("💧 遗失的魔术贴 (Missing Velcro)")
    st.markdown("挖掘 **高 IDR (无序区)** 驱动的、可能被传统结构预测遗漏的相分离互作。")
    
    tab1, tab2 = st.tabs(["🔍 单蛋白分析 (需 metapredict)", "💎 全局挖掘 (基于数据)"])
    
    with tab1:
        g_in = st.text_input("输入基因名:", "FUS")
        if st.button("计算序列 IDR"):
            try:
                import metapredict as meta
                import requests
                # 简化版逻辑
                url = "https://rest.uniprot.org/uniprotkb/search"
                res = requests.get(url, params={'query': f"gene_exact:{g_in} AND reviewed:true", 'format': 'json'}).json()
                if res['results']:
                    seq = res['results'][0]['sequence']['value']
                    scores = meta.predict_disorder(seq)
                    frac = sum(1 for s in scores if s>0.5)/len(seq)
                    st.metric("IDR 比例", f"{frac:.2%}")
                    st.line_chart(scores)
                else:
                    st.error("未找到序列")
            except:
                st.warning("云端缺少 metapredict 或网络问题。")
                
    with tab2:
        if st.button("扫描潜在液滴互作"):
            # 扫描图谱中 IDR > 0.4 的节点
            high_idr_nodes = []
            for n, d in G.nodes(data=True):
                # 尝试从节点属性读取 idr (如果之前存过)
                idr = d.get('idr', -1)
                if idr > 0.4:
                    high_idr_nodes.append({'Gene': n, 'IDR': idr, 'Loc': d.get('loc')})
            
            if not high_idr_nodes:
                st.warning("当前图谱数据中未包含预计算的 IDR 信息。请在 Colab 中运行 '批量计算' 并保存到 CSV。")
            else:
                st.write(f"发现 {len(high_idr_nodes)} 个高无序蛋白。正在两两配对...")
                # 简单的两两配对逻辑 (取前 50 个演示)
                cands = sorted(high_idr_nodes, key=lambda x:x['IDR'], reverse=True)[:50]
                pairs = []
                for i in range(len(cands)):
                    for j in range(i+1, len(cands)):
                        u, v = cands[i], cands[j]
                        # 如果没有已知互作，且 AI 分数低，但都有高 IDR -> 可能是 LLPS
                        if not G.has_edge(u['Gene'], v['Gene']):
                            pairs.append({
                                "Protein A": u['Gene'], "IDR A": u['IDR'],
                                "Protein B": v['Gene'], "IDR B": v['IDR'],
                                "Prediction": "💧 LLPS Potential"
                            })
                
                if pairs:
                    st.dataframe(pd.DataFrame(pairs))
                else:
                    st.info("未发现显著配对。")
