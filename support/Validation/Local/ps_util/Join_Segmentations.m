
Segment_bpt_MS=Sbpt_MS(:,:,9);

labels = unique(Segment_bpt_MS);
 Segment_new_MS=zeros(size(Skmeans_MS));
 for kk=1:size(Skmeans_MS,3),
 Segment_kmeans_MS=Skmeans_MS(:,:,kk);
  Segment_new_kk=zeros(size(Segment_new_MS,1),size(Segment_new_MS,2));
for kkk=1:length(labels)
   
    idx = Segment_bpt_MS==labels(kkk);
    label_vector=Segment_kmeans_MS(idx);
%     [a,b] = histc(label_vector,unique(label_vector));
    [a,b] = hist(label_vector,unique(label_vector));
%     y=a(b);
    [max_a,index_a] = max(a);
    Segment_new_kk(idx)=b(index_a);
end
Segment_new_MS(:,:,kk)=Segment_new_kk;
 end
 
 Segment_bpt_PAN=Sbpt_MS(:,:,9);

labels = unique(Segment_bpt_PAN);
 Segment_new_PAN=zeros(size(Skmeans_PAN));
 for kk=1:size(Skmeans_PAN,3),
 Segment_kmeans_PAN=Skmeans_PAN(:,:,kk);
  Segment_new_kk=zeros(size(Segment_new_PAN,1),size(Segment_new_PAN,2));
for kkk=1:length(labels)
   
    idx = Segment_bpt_PAN==labels(kkk);
    label_vector=Segment_kmeans_PAN(idx);
%     [a,b] = histc(label_vector,unique(label_vector));
    [a,b] = hist(label_vector,unique(label_vector));
%     y=a(b);
    [max_a,index_a] = max(a);
    Segment_new_kk(idx)=b(index_a);
end
Segment_new_PAN(:,:,kk)=Segment_new_kk;
 end