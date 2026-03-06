TH2D* SwapAxesManually(TH2D* h_original) {
    if (!h_original) return nullptr;
    
    // 获取原始直方图的参数
    int nx = h_original->GetNbinsX();
    int ny = h_original->GetNbinsY();
    double xmin = h_original->GetXaxis()->GetXmin();
    double xmax = h_original->GetXaxis()->GetXmax();
    double ymin = h_original->GetYaxis()->GetXmin();
    double ymax = h_original->GetYaxis()->GetXmax();
    
    // 获取标题
    TString title = h_original->GetTitle();
    TString xtitle = h_original->GetXaxis()->GetTitle();
    TString ytitle = h_original->GetYaxis()->GetTitle();
    
    // 创建新直方图（交换轴）
    TH2D* h_swapped = new TH2D(
        Form("%s_swapped", h_original->GetName()),
        Form("%s (swapped)", title.Data()),
        ny, ymin, ymax,  // 原Y轴变为新X轴
        nx, xmin, xmax   // 原X轴变为新Y轴
    );
    
    // 交换轴标题
    h_swapped->GetXaxis()->SetTitle(ytitle);
    h_swapped->GetYaxis()->SetTitle(xtitle);
    
    // 复制数据（转置）
    for (int ix = 1; ix <= nx; ix++) {
        for (int iy = 1; iy <= ny; iy++) {
            double content = h_original->GetBinContent(ix, iy);
            double error = h_original->GetBinError(ix, iy);
            
            // 转置：(ix, iy) -> (iy, ix)
            h_swapped->SetBinContent(iy, ix, content);
            h_swapped->SetBinError(iy, ix, error);
        }
    }
    
    // 复制统计信息
    h_swapped->SetEntries(h_original->GetEntries());
    
    return h_swapped;
}

