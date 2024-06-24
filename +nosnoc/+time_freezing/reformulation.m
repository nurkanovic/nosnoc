function model = reformulation(cls_model, opts)
    if any(cls_model.e ~= 0)
        model = nosnoc.time_freezing.cls_elastic(cls_model, opts);
    else
        if opts.tf_multicontact
            model = nosnoc.time_freezing.cls_inelastic_multicontact(cls_model, opts);
        else
            model = nosnoc.time_freezing.cls_inelastic(cls_model, opts);
        end
    end
end

